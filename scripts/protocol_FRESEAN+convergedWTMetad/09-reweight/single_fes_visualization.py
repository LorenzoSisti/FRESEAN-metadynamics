#!/usr/bin/env python3
"""
single_fes_visualization.py
----------------------------
Genera una Free Energy Surface (FES) 2D a partire dal reweighting di una
metadinamica, combinando:
  - plumed-output.out : time, pca.eig-1, pca.eig-2, pca.residual, metad.bias
  - gyrate.xvg         : raggio di girazione (Rg) da gmx gyrate
  - sasa.xvg            : area SASA totale da gmx sasa

Lo script:
  1. Carica e allinea i tre file (tutti sullo stesso set di frame).
  2. Salva un unico file aggregato (utile per riusi futuri / altri plot).
  3. Calcola i pesi di reweighting w_i = exp(bias_i / kT)  (Tiwary-Parrinello,
     valido quando il bias e' quasi-statico, es. fine di una WT-metadynamics).
  4. Costruisce un istogramma 2D pesato sulle CV scelte (default: pca.eig-1
     vs pca.eig-2, per riprodurre il tipo di mappa mostrato nell'esempio
     KRAS) e lo converte in energia libera G = -kT ln(P).
  5. Plotta la FES con matplotlib in stile "jet" con cap a ENERGY_CAP e
     grigio per le regioni non campionate / oltre soglia.

Posizione attesa dei file (secondo la tua struttura):
  09-reweight/single_fes_visualization.py   <- questo script
  09-reweight/reweight_0/plumed-output.out
  09-reweight/reweight_0/gyrate.xvg
  09-reweight/reweight_0/sasa.xvg

Uso:
  cd 09-reweight
  python single_fes_visualization.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ======================================================================
# CONFIGURAZIONE
# ======================================================================
DATA_DIR = "reweight_0"
PLUMED_FILE = os.path.join(DATA_DIR, "plumed-output.out")
GYRATE_FILE = os.path.join(DATA_DIR, "gyrate.xvg")
SASA_FILE = os.path.join(DATA_DIR, "sasa.xvg")

OUTPUT_DAT = "reweight_aggregated.dat"     # file aggregato con tutte le colonne
OUTPUT_FIG = "FES_pca1_pca2.png"           # figura finale

# Costante di Boltzmann e temperatura della simulazione (adatta se necessario)
kB = 0.0083144621   # kJ/(mol*K)
T = 300.0           # K
kBT = kB * T

# Quali colonne usare come assi della FES:
#   "pca1", "pca2", "residual", "rg", "sasa"
CV_X = "rg"
CV_Y = "sasa"

NBINS = 100          # numero di bin per asse dell'istogramma 2D

# ======================================================================
# 1. CARICAMENTO DATI
# ======================================================================
print("Caricamento dei file...")

plumed_data = np.loadtxt(PLUMED_FILE, skiprows=1)
time_plumed = plumed_data[:, 0]
pca1 = plumed_data[:, 1]
pca2 = plumed_data[:, 2]
residual = plumed_data[:, 3]
bias = plumed_data[:, 4]  # bias di metadinamica per il reweighting

gyrate_data = np.loadtxt(GYRATE_FILE, comments=("#", "@"))
rg = gyrate_data[:, 1]  # Rg totale (colonna 1)

sasa_data = np.loadtxt(SASA_FILE, comments=("#", "@"))
sasa = sasa_data[:, 1]  # SASA totale

# ------------------------------------------------------------------
# Allineamento: tronca tutti gli array alla lunghezza minima comune,
# nel caso in cui gmx e plumed abbiano scritto un numero leggermente
# diverso di frame (es. per l'ultimo frame non completato).
# ------------------------------------------------------------------
n = min(len(bias), len(rg), len(sasa))
if not (len(bias) == len(rg) == len(sasa)):
    print(f"Attenzione: lunghezze diverse "
          f"(plumed={len(bias)}, gyrate={len(rg)}, sasa={len(sasa)}). "
          f"Taglio tutto a {n} frame.")
else:
    print("Tutti i file sono caricati e allineati correttamente!")

time_plumed = time_plumed[:n]
pca1 = pca1[:n]
pca2 = pca2[:n]
residual = residual[:n]
bias = bias[:n]
rg = rg[:n]
sasa = sasa[:n]

# ======================================================================
# 2. SALVATAGGIO FILE AGGREGATO
# ======================================================================
header = "time(ps)  pca.eig-1  pca.eig-2  pca.residual  metad.bias(kJ/mol)  Rg(nm)  SASA(nm^2)"
np.savetxt(
    OUTPUT_DAT,
    np.column_stack([time_plumed, pca1, pca2, residual, bias, rg, sasa]),
    header=header,
    fmt="%.6f",
)
print(f"File aggregato salvato in: {OUTPUT_DAT}")

# ======================================================================
# 3. PESI DI REWEIGHTING
# ======================================================================
# w_i = exp(V(s_i,t) / kT)  -- valido per bias quasi-statico (fine WT-MetaD)
weights = np.exp(bias / kBT)
weights /= weights.sum()

# ======================================================================
# 4. ISTOGRAMMA 2D PESATO -> FES
# ======================================================================
cv_map = {
    "pca1": (pca1, "PC1"),
    "pca2": (pca2, "PC2"),
    "residual": (residual, "PCA residual"),
    "rg": (rg, "Rg (nm)"),
    "sasa": (sasa, r"SASA (nm$^2$)"),
}
cv_x, label_x = cv_map[CV_X]
cv_y, label_y = cv_map[CV_Y]

hist, xedges, yedges = np.histogram2d(cv_x, cv_y, bins=NBINS, weights=weights)

# Evita log(0): i bin mai visitati diventano NaN (mostrati in grigio)
hist_masked = np.where(hist > 0, hist, np.nan)

# G(x,y) = -kT ln[P(x,y)], riportata al minimo (G_min = 0)
fes = -kBT * np.log(hist_masked)
fes -= np.nanmin(fes)
fes = fes.T  # orientamento corretto per pcolormesh (righe=y, colonne=x)

ENERGY_CAP = np.ceil(np.nanmax(fes))

x_centers = 0.5 * (xedges[:-1] + xedges[1:])
y_centers = 0.5 * (yedges[:-1] + yedges[1:])

# ======================================================================
# 5. PLOT
# ======================================================================
cmap = plt.cm.jet.copy()
cmap.set_bad(color="lightgray")    # bin non campionati (NaN)
cmap.set_over(color="lightgray")   # energie oltre ENERGY_CAP

norm = mcolors.Normalize(vmin=0, vmax=ENERGY_CAP)

fig, ax = plt.subplots(figsize=(6, 5))

pcm = ax.pcolormesh(x_centers, y_centers, fes, cmap=cmap, norm=norm, shading="gouraud")

cbar = fig.colorbar(pcm, ax=ax, extend="max")
cbar.set_label(r"$\Delta G$ (kJ/mol)", fontsize=12)

ax.set_xlabel(label_x, fontsize=12)
ax.set_ylabel(label_y, fontsize=12)
ax.set_title("FES", fontsize=14)

plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=300)
print(f"Figura salvata in: {OUTPUT_FIG}")

plt.show()

print("Max FES:", np.nanmax(fes))
print("99° percentile:", np.nanpercentile(fes, 99))
print("95° percentile:", np.nanpercentile(fes, 95))
