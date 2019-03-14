import cffi
import gzip


ffi = cffi.FFI()
ffi.cdef("""
    typedef void* variant;

    variant variant_new(
        char *name,
        char *chrom,
        unsigned int position,
        char *a1,
        char *a2
    );
    void* variant_free(variant);
    void* variant_print(variant);
    void* variant_complement_alleles(variant);
""")

# C = ffi.dlopen("../target/debug/librsgeneparselib.dylib")
C = ffi.dlopen("../target/release/librsgeneparselib.dylib")


class Variant(object):
    __slots__ = ["_obj"]

    def __init__(self, name, chrom, position, alleles):
        if type(name) is bytes:
            self.init_from_bytes(name, chrom, position, alleles)

        else:
            self.init_from_bytes(name, chrom, position, alleles)

    def init_from_str(self, name, chrom, position, alleles):
        a1, a2 = alleles

        name = ffi.new("char[]", name.encode())
        chrom = ffi.new("char[]", str(chrom).encode())
        a1 = ffi.new("char[]", a1.encode())
        a2 = ffi.new("char[]", a2.encode())

        self._obj = C.variant_new(name, chrom, position, a1, a2)

    def init_from_bytes(self, name, chrom, position, alleles):
        a1, a2 = alleles
        self._obj = C.variant_new(name, chrom, position, a1, a2)

    def complement_alleles(self):
        C.variant_complement_alleles(self._obj)

    def print(self):
        C.variant_print(self._obj)

    def __del__(self):
        C.variant_free(self._obj)


if __name__ == "__main__":
    filename = "/Users/legaultmarc/projects/StatGen/lvef/full_lvef_uk_biobank_gwas_results.txt.gz"

    variants = []
    with gzip.open(filename, "rb") as f:
        header = next(f).decode("utf-8").strip().split()
        header = {col: idx for idx, col in enumerate(header)}

        name_col = header.get("snp")

        for line in f:
            line = line.strip().split()


            v = Variant(
                line[name_col] if name_col else b"",
                line[header["chr"]],
                int(line[header["pos"]]),
                [line[header["major"]],
                 line[header["minor"]]],
            )

            variants.append(v)
