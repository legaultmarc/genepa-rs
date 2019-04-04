🚧🚧🚧🚧🚧🚧🚧

**This is currently experimental / for fun / a work in progress**

🚧🚧🚧🚧🚧🚧🚧

# Background

We wrote [geneparse](https://github.com/pgxcentre/geneparse) which provides a very convenient API to represent variants
and genotypes in a unified format. We are also adding support for different
genotype file formats which makes it easy to develop tools that are very
versatile.

The library is currently more or less pure Python which has two major
drawbacks:

1. It is relatively slow
2. It is hard to use in other languages (_e.g._ in R or C/C++)

The current project is mostly an experiment, but it is a reimplementation of
`geneparse` in Rust. Rust has the features of modern programming languages and
it is very robust and fast. We aim to gradually implement the "workhorse"
classes and parsers from geneparse and to use them as drop-in replacements
from the pure Python implementations.

# How does it work

Basically, there are three major components:

1. The pure Rust implementations that could be used natively (core.rs)
2. A C compatibility layer (c_api.rs)
3. The Python side which uses CFFI to call the C API (python/api.py)

To use this as a drop in replacement would be very easy. The only necessary
work would be to add ``__getattribute__`` and ``__setattribute__`` methods
on the Python side and bind them to getters and setters on the Rust API.
Then, simply inherting from ``geneparse.Variant`` would take care of the currently unimplemented methods.

# Benchmarks

Now that there is a draft Plink reader implementation, I have done some
benchmarking. I read a Plink binary file and printed the variants as well as
the computed MAF.

The BED file was 2.7GB containing 1,705,969 variants and 6,399 individuals.

The Rust implementation took 6m42.094s for this task when creating the bim
index and 6m31.473s when repeating the task after creating the index.

The Python implementation took 12m55.804s (an almost 2x speedup).

# How to run it

For now, this is very hackish. You need to manually compile the library using
``cargo build`` or ``cargo build --release``. This will create the dylib, so
or dll that is required for CFFI to find the relevant functions. The path
is currently hard coded so it will not work on platforms other than MacOS.
Then from the ``python/`` directory you can run the ``api.py`` script.

# Acknowledgements

I used this blog post to better understand most of the FFI machinery that I
use in here:

[https://bheisler.github.io/post/calling-rust-in-python/](https://bheisler.github.io/post/calling-rust-in-python/)

