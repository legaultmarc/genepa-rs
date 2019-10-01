import cffi
import gzip

import numpy as np


ffi = cffi.FFI()
ffi.cdef("""
    typedef void* variant;
    typedef void* genotypes;
    typedef void* plink_reader;

    // Variant
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

    // PlinkReader
    plink_reader plink_reader_new(char *prefix);
    void* plink_reader_free(plink_reader);
    genotypes plink_reader_next(plink_reader);

    // Genotypes
    void* genotypes_print(genotypes);
    variant genotypes_get_variant(genotypes);
    float* genotypes_get_genotypes(genotypes);
    unsigned int genotypes_len(genotypes);
    double genotypes_maf(genotypes);
    void* genotypes_free(genotypes);
    void* genotypes_data_free(genotypes);

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

    @classmethod
    def new_from_pointer(cls, obj):
        x = cls.__new__(cls)
        x._obj = obj
        return x

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


class Genotypes(object):
    __slots__ = ["_obj", "_data_pointer"]

    def __del__(self):
        C.genotypes_free(self._obj)

        if self._data_pointer is not None:
            C.genotypes_data_free(self._data_pointer)

    @classmethod
    def new_from_pointer(cls, obj):
        x = cls()
        x._obj = obj
        x._data_pointer = None
        return x

    @property
    def variant(self):
        return Variant.new_from_pointer(C.genotypes_get_variant(self._obj))

    @property
    def genotypes(self):
        n = len(self)

        if self._data_pointer is None:
            self._data_pointer = C.genotypes_get_genotypes(self._obj)

        return np.frombuffer(
            ffi.buffer(self._data_pointer, 4 * n),
            count=n,
            dtype=np.float32
        )

    def print(self):
        C.genotypes_print(self._obj)

    def maf(self):
        return C.genotypes_maf(self._obj)

    def __len__(self):
        return C.genotypes_len(self._obj)


def print_genotypes(obj):
    C.genotypes_print(obj)


class PlinkReader(object):
    __slots__ = ["_obj"]

    def __init__(self, prefix):
        self._obj = C.plink_reader_new(prefix.encode("ascii"))

    def __del__(self):
        C.plink_reader_free(self._obj)

    def __iter__(self):
        return self

    def __next__(self):
        ptr = C.plink_reader_next(self._obj)

        if ptr == ffi.NULL:
            raise StopIteration()

        return Genotypes.new_from_pointer(ptr)


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
