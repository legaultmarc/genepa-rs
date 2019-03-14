import gzip

from geneparse import Variant

if __name__ == "__main__":
    filename = "/Users/legaultmarc/projects/StatGen/lvef/full_lvef_uk_biobank_gwas_results.txt.gz"

    variants = []
    with gzip.open(filename, "rt") as f:
        header = {col: idx for idx, col in enumerate(next(f).strip().split())}

        name_col = header.get("snp")

        for line in f:
            line = line.strip().split()


            v = Variant(
                line[name_col] if name_col else "",
                line[header["chr"]],
                int(line[header["pos"]]),
                [line[header["major"]],
                 line[header["minor"]]],
            )

            variants.append(v)
