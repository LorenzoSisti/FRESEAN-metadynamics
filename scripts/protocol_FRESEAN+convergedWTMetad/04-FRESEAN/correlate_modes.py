#!/usr/bin/env python3
"""
correlate_modes.py

Calcola la correlazione di Pearson fra i modi vibrazionali a frequenza
zero (tipicamente modi 7, 8, 9) estratti indipendentemente da repliche
diverse della stessa simulazione, replicando la logica della Fig. 1 di
Sauer et al., Sci. Adv. 12, eaea4617 (2026).

Nessuna dipendenza esterna: solo libreria standard Python (math, sys,
itertools, argparse). Non richiede numpy.

INPUT ATTESO
------------
File .xyz multi-frame prodotti da `fresean extract` (es.
evec_freq1_mode1-30_cg.xyz), con formato standard XYZ ripetuto:

    <numero atomi>
    <riga di commento>
    X  dx1  dy1  dz1
    X  dx2  dy2  dz2
    ...
    <numero atomi>
    <riga di commento>
    X  dx1  dy1  dz1
    ...

dove ogni blocco ("frame") corrisponde a un modo, in ordine crescente
a partire da modeStart (di norma 1, come da extract.inp del protocollo).

NOTA SUL SEGNO DEGLI AUTOVETTORI
---------------------------------
Un autovettore e il suo opposto (-v) rappresentano la stessa direzione
fisica di moto collettivo: la diagonalizzazione non fissa il segno.
Per questo motivo, come nel paper, viene riportato il VALORE ASSOLUTO
della correlazione di Pearson fra due modi, non il valore con segno.

USO
---
python3 correlate_modes.py \
    --modes 7 8 9 \
    R1/evec_freq1_mode1-30_cg.xyz \
    R2/evec_freq1_mode1-30_cg.xyz \
    R3/evec_freq1_mode1-30_cg.xyz \
    R4/evec_freq1_mode1-30_cg.xyz \
    R5/evec_freq1_mode1-30_cg.xyz

Il primo file passato è preso come riferimento (analogo a R1 nel paper):
vengono stampate le matrici di correlazione fra il riferimento e ciascuno
degli altri file, per tutte le combinazioni fra i modi richiesti.
"""

import sys
import argparse
import math
from itertools import combinations


def parse_xyz_modes(filename, mode_start=1):
    """
    Legge un file .xyz multi-frame e restituisce un dizionario
    {indice_modo: [dx1, dy1, dz1, dx2, dy2, dz2, ...]}

    indice_modo parte da mode_start (default 1), coerente con il
    parametro #modeStart usato in extract.inp dal protocollo FRESEAN.
    """
    modes = {}
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip() != ""]

    i = 0
    mode_idx = mode_start
    n_lines = len(lines)
    while i < n_lines:
        try:
            natoms = int(lines[i])
        except ValueError:
            raise ValueError(
                f"Riga {i+1} di {filename} non e' un intero valido "
                f"(atteso il numero di atomi di un frame): '{lines[i]}'"
            )
        i += 1  # salta la riga con natoms
        i += 1  # salta la riga di commento

        coords = []
        for _ in range(natoms):
            if i >= n_lines:
                raise ValueError(
                    f"File {filename} troncato: attesi {natoms} atomi "
                    f"per il modo {mode_idx}, trovati meno."
                )
            parts = lines[i].split()
            if len(parts) < 4:
                raise ValueError(
                    f"Riga {i+1} di {filename} malformata, attese almeno "
                    f"4 colonne (elemento x y z): '{lines[i]}'"
                )
            # parts[0] e' l'etichetta (tipicamente "X"), parts[1:4] sono x,y,z
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            coords.extend([x, y, z])
            i += 1

        modes[mode_idx] = coords
        mode_idx += 1

    return modes


def check_translation_character(mode_vector, natoms):
    """
    Verifica se un modo ha carattere traslazionale puro: in tal caso
    ogni atomo si sposta nella stessa identica direzione, cioe' il
    vettore spostamento di ciascun atomo e' parallelo al vettore
    spostamento medio dell'intero modo (cosine similarity vicina a 1).

    Restituisce (cos_medio, cos_minimo) sui singoli atomi rispetto al
    vettore medio. Valori entrambi vicini a 1.0 (es. > 0.95) indicano
    un modo traslazionale; valori bassi o dispersi indicano che NON e'
    un modo di traslazione di corpo rigido (atteso per i modi 7+).
    """
    per_atom = [mode_vector[3*i:3*i+3] for i in range(natoms)]
    mean_vec = [
        sum(a[k] for a in per_atom) / natoms for k in range(3)
    ]
    mean_norm = math.sqrt(sum(x*x for x in mean_vec))
    if mean_norm == 0.0:
        return 0.0, 0.0

    coss = []
    for a in per_atom:
        a_norm = math.sqrt(sum(x*x for x in a))
        if a_norm == 0.0:
            coss.append(0.0)
            continue
        dot = sum(a[k]*mean_vec[k] for k in range(3))
        coss.append(dot / (a_norm * mean_norm))

    return sum(coss)/len(coss), min(coss)


def sanity_check_modes_1_to_6(filename, natoms_hint=None):
    """
    Legge i modi 1-6 (attesi: traslazione 1-3, rotazione 4-6) da un file
    .xyz e stampa un controllo diagnostico. Da eseguire PRIMA di fidarsi
    che i modi 7/8/9 siano davvero le vibrazioni anarmoniche a maggior
    contributo, escluse traslazione/rotazione, come atteso dal protocollo.
    """
    modes = parse_xyz_modes(filename, mode_start=1)
    natoms = len(modes[1]) // 3
    print(f"\n=== Sanity check modi 1-6 (traslazione/rotazione attese) su {filename} ===")
    print(f"Numero di atomi (coarse-grained) rilevato: {natoms}")
    for m in range(1, 7):
        if m not in modes:
            print(f"  Modo {m}: NON TROVATO nel file (file troppo corto?)")
            continue
        avg_cos, min_cos = check_translation_character(modes[m], natoms)
        label = "traslazione pura (atteso 1-3)" if avg_cos > 0.95 else \
                "NON traslazionale (atteso per rotazione 4-6, o anomalia se e' modo 1-3)"
        print(f"  Modo {m}: cos medio con vettore medio = {avg_cos:.3f}, "
              f"cos minimo = {min_cos:.3f}  -> {label}")
    print("Se i modi 1-3 non risultano chiaramente traslazionali (cos medio "
          "vicino a 1), la numerazione/ordinamento dei modi va verificata "
          "prima di procedere con l'analisi dei modi 7/8/9.\n")


def cosine_similarity(v1, v2):
    """
    Similarita' del coseno pura Python fra due vettori (liste di float)
    della stessa lunghezza: prodotto scalare normalizzato, in [-1, 1].

    Usata al posto della correlazione di Pearson perche' i modi FRESEAN
    sono vettori di spostamento fisico: sottrarre la media delle componenti
    (come fa Pearson) non ha significato fisico qui. Il coseno fra vettori
    normalizzati e' la metrica corretta per confrontare due autovettori.
    """
    n = len(v1)
    if n != len(v2):
        raise ValueError(
            f"I due modi hanno un numero diverso di componenti "
            f"({n} vs {len(v2)}): probabile mismatch nel numero di "
            f"atomi coarse-grained fra le repliche confrontate."
        )
    if n == 0:
        raise ValueError("Vettori vuoti, impossibile calcolare la similarita'.")

    dot = sum(a * b for a, b in zip(v1, v2))
    norm1 = math.sqrt(sum(a * a for a in v1))
    norm2 = math.sqrt(sum(b * b for b in v2))

    if norm1 == 0.0 or norm2 == 0.0:
        sys.stderr.write(
            "ATTENZIONE: uno dei due modi confrontati ha norma nulla; "
            "similarita' impostata a 0.0 (verificare il file di input).\n"
        )
        return 0.0

    return dot / (norm1 * norm2)


def normalize(v):
    """Restituisce v scalato a norma unitaria (lista di float)."""
    norm = math.sqrt(sum(a * a for a in v))
    if norm == 0.0:
        raise ValueError("Impossibile normalizzare un vettore di norma nulla.")
    return [a / norm for a in v]


def determinant(matrix):
    """
    Determinante di una matrice quadrata (lista di liste), via espansione
    di Laplace. Adatto solo a matrici piccole (qui: 2x2 o 3x3), che e'
    esattamente il caso d'uso (sottospazi a 2 o 3 modi). Puro Python,
    nessuna dipendenza esterna.
    """
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    det = 0.0
    for col in range(n):
        minor = [row[:col] + row[col+1:] for row in matrix[1:]]
        sign = 1 if col % 2 == 0 else -1
        det += sign * matrix[0][col] * determinant(minor)
    return det


def subspace_correlation(ref_modes, other_modes, mode_indices):
    """
    Correlazione fra il sottospazio spannato da un insieme di modi nella
    replica di riferimento e il sottospazio spannato dagli stessi indici
    di modo nell'altra replica (analogo alle correlazioni 2D/3D di Fig. 1,
    eq. S13/S14 del paper). Necessaria perche' confrontare solo modo-con-
    modo (diagonale) sottostima la riproducibilita' quando due modi quasi
    degeneri (es. 8 e 9) si scambiano di ordine fra repliche diverse.

    Calcolata come |det(C)|, dove C_ij = coseno fra il modo mode_indices[i]
    del riferimento e il modo mode_indices[j] dell'altra replica, con tutti
    i vettori pre-normalizzati. Se i due insiemi di modi fossero identici
    (stessa base ortonormale), |det(C)| = 1; quanto piu' i sottospazi
    differiscono, tanto piu' il valore scende verso 0.
    """
    n = len(mode_indices)
    C = [[0.0] * n for _ in range(n)]
    for i, mi in enumerate(mode_indices):
        u = normalize(ref_modes[mi])
        for j, mj in enumerate(mode_indices):
            v = normalize(other_modes[mj])
            C[i][j] = sum(a * b for a, b in zip(u, v))
    return abs(determinant(C))


def print_matrix(ref_label, other_label, mode_list, ref_modes, other_modes):
    """Stampa una matrice di correlazione |Pearson| fra mode_list x mode_list,
    fra il file di riferimento (ref) e un altro file (other)."""
    header = "        " + "".join(f"{m:>10}" for m in mode_list)
    print(f"\n=== {ref_label} (righe) vs {other_label} (colonne) ===")
    print(header)
    for m_ref in mode_list:
        row_vals = []
        for m_other in mode_list:
            if m_ref not in ref_modes:
                raise ValueError(f"Modo {m_ref} non trovato in {ref_label}")
            if m_other not in other_modes:
                raise ValueError(f"Modo {m_other} non trovato in {other_label}")
            r = cosine_similarity(ref_modes[m_ref], other_modes[m_other])
            row_vals.append(abs(r))
        row_str = "".join(f"{v:10.3f}" for v in row_vals)
        print(f"mode {m_ref:>3}{row_str}")


def main():
    parser = argparse.ArgumentParser(
        description="Correlazione (coseno + sottospazio) fra modi FRESEAN di repliche "
                    "diverse (analogo alla Fig. 1 e Fig. S5 di Sauer et al., Sci. Adv. 2026)."
    )
    parser.add_argument(
        "files", nargs="+",
        help="File .xyz multi-frame (uno per replica)."
    )
    parser.add_argument(
        "--modes", type=str, default="7,8,9",
        help="Indici dei modi da confrontare, separati da virgola, "
             "senza spazi (default: 7,8,9). Esempio: --modes 7,8,9"
    )
    parser.add_argument(
        "--mode-start", type=int, default=1,
        help="Indice del primo frame nel file .xyz (default: 1, coerente con "
             "modeStart=1 in extract.inp)."
    )
    parser.add_argument(
        "--skip-sanity-check", action="store_true",
        help="Salta il controllo diagnostico sui modi 1-6 (traslazione/rotazione)."
    )
    args = parser.parse_args()

    if len(args.files) < 2:
        print("Servono almeno 2 file da confrontare.")
        sys.exit(1)

    try:
        mode_list = [int(x) for x in args.modes.split(",")]
    except ValueError:
        print(f"Errore: --modes deve essere una lista di interi separati da "
              f"virgola senza spazi, es. 7,8,9. Ricevuto: '{args.modes}'")
        sys.exit(1)
    args.modes = mode_list

    print(f"Modi confrontati: {args.modes}")

    parsed = {}
    for fname in args.files:
        print(f"Lettura {fname} ...")
        parsed[fname] = parse_xyz_modes(fname, mode_start=args.mode_start)

    # --- Sanity check preliminare: i modi 1-6 sono davvero trasl./rot.? ---
    if not args.skip_sanity_check:
        if args.mode_start > 1:
            print("\n(--mode-start > 1: salto il sanity check sui modi 1-6, "
                  "non presenti nei file forniti)")
        else:
            for fname in args.files:
                sanity_check_modes_1_to_6(fname)

    # --- Correlazione su TUTTE le coppie di repliche, non solo la prima ---
    pairs = list(combinations(args.files, 2))
    pair_2d_results = {}   # (fileA, fileB) -> valore subspace 2D
    pair_nd_results = {}   # (fileA, fileB) -> valore subspace ND (tutti i modi richiesti)

    for fA, fB in pairs:
        modesA = parsed[fA]
        modesB = parsed[fB]
        print_matrix(fA, fB, args.modes, modesA, modesB)

        if len(args.modes) >= 2:
            pair_2d_results[(fA, fB)] = subspace_correlation(modesA, modesB, args.modes[:2])
        if len(args.modes) >= 3:
            pair_nd_results[(fA, fB)] = subspace_correlation(modesA, modesB, args.modes)

    # --- Riepilogo: correlazione di sottospazio per TUTTE le coppie ---
    if pair_2d_results:
        print(f"\n=== Correlazione di sottospazio 2D (modi {args.modes[:2]}), tutte le coppie ===")
        for (fA, fB), val in pair_2d_results.items():
            print(f"{fA}  vs  {fB}:  |det| = {val:.3f}")

    if pair_nd_results:
        print(f"\n=== Correlazione di sottospazio {len(args.modes)}D (modi {args.modes}), tutte le coppie ===")
        for (fA, fB), val in pair_nd_results.items():
            print(f"{fA}  vs  {fB}:  |det| = {val:.3f}")

    # --- Classifica: per ciascuna replica, media della sua correlazione
    #     di sottospazio ND con TUTTE le altre 4 -> guida la scelta del
    #     riferimento (non piu' assunto a priori come il primo file) ---
    results_for_ranking = pair_nd_results if pair_nd_results else pair_2d_results
    if results_for_ranking:
        dim_label = f"{len(args.modes)}D" if pair_nd_results else "2D"
        per_file_scores = {f: [] for f in args.files}
        for (fA, fB), val in results_for_ranking.items():
            per_file_scores[fA].append(val)
            per_file_scores[fB].append(val)

        print(f"\n=== Classifica repliche per correlazione media di sottospazio {dim_label} "
              f"(con TUTTE le altre repliche) ===")
        ranking = sorted(
            per_file_scores.items(),
            key=lambda kv: sum(kv[1])/len(kv[1]),
            reverse=True
        )
        for fname, scores in ranking:
            avg = sum(scores)/len(scores)
            print(f"{fname}:  media = {avg:.3f}   (valori: "
                  + ", ".join(f"{v:.3f}" for v in scores) + ")")

        best = ranking[0][0]
        print(f"\nReplica suggerita come riferimento (correlazione media piu' alta "
              f"con le altre 4): {best}")


if __name__ == "__main__":
    main()
