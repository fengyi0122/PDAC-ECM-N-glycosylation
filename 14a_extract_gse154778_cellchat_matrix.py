import argparse
import gzip
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--output-matrix", required=True)
    parser.add_argument("--output-library", required=True)
    parser.add_argument("--output-audit", required=True)
    args = parser.parse_args()

    wanted = []
    with open(args.genes, encoding="utf-8") as handle:
        handle.readline()
        for line in handle:
            gene = line.rstrip("\r\n").split("\t")[0]
            if gene:
                wanted.append(gene)
    wanted_set = set(wanted)

    selected = {}
    with gzip.open(args.input, "rt", encoding="utf-8-sig", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        cells = header[1:]
        library = np.zeros(len(cells), dtype=np.float64)
        for line in handle:
            split_at = line.find(",")
            if split_at < 1:
                continue
            gene = line[:split_at]
            values = np.fromstring(line[split_at + 1 :], sep=",", dtype=np.float64)
            if values.size != len(cells):
                raise RuntimeError(f"Column count mismatch for {gene}: {values.size} versus {len(cells)}")
            library += values
            if gene in wanted_set:
                if gene in selected:
                    selected[gene] += values
                else:
                    selected[gene] = values.copy()

    Path(args.output_matrix).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output_matrix, "wt", encoding="utf-8", newline="") as handle:
        handle.write("gene\t" + "\t".join(cells) + "\n")
        for gene in wanted:
            values = selected.get(gene, np.zeros(len(cells), dtype=np.float64))
            handle.write(gene + "\t" + "\t".join(np.char.mod("%.10g", values)) + "\n")

    with open(args.output_library, "w", encoding="utf-8", newline="") as handle:
        handle.write("cell_id\tlibrary_size\n")
        for cell, value in zip(cells, library):
            handle.write(f"{cell}\t{value:.10g}\n")

    with open(args.output_audit, "w", encoding="utf-8", newline="") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"cells\t{len(cells)}\n")
        handle.write(f"requested_genes\t{len(wanted)}\n")
        handle.write(f"recovered_genes\t{len(selected)}\n")
        handle.write(f"positive_library_cells\t{int(np.sum(library > 0))}\n")


if __name__ == "__main__":
    main()
