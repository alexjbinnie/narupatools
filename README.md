# narupatools Python Package

This package includes extensions and utilities for working with the Narupa framework.

# Installation

This library can be installed using

```conda install narupatools -c conda-forge -c omnia -c irl -c alexjbinnie```

# Disclaimer

narupatools is a separate project which is not a part of Narupa, and does not have any affiliation with either Narupa or
the Intangible Realities Laboratory. The author would like to make it clear that issues arising due to this library
should not be reported on the Narupa repositories, and instead reported on this repository.

# License

This project is released under the GNU General Public License version 3. All code written for narupatools is (c)
Alex Jamieson-Binnie, All Rights Reserved, unless otherwise specified in the source file.

# Development

## Testing

Unit tests for narupatools can be found in `narupatools/tests`, and can be run
using [pytest](https://docs.pytest.org/en/stable/). The [Hypothesis](https://hypothesis.readthedocs.io/en/latest/)
library is also used for testing, to generate test data for certain methods.

## Static Analysis

All code in this library has type annotations, which are checked by [mypy](https://github.com/python/mypy). Type
annotations are strictly enforced.

Stubs for packages referenced by this project that do not have their own type annotations can be found in the `stubs/`
directory. This currently contains type annotations for the methods, classes and variables used by `narupatools` in ASE,
matplotlib, MDAnalysis, mdtraj, nglview and OpenMM.

The [flake8](https://flake8.pycqa.org/en/latest/) tool is also run on the source, tests and stubs to enforce a set of
style guidelines. In addition to the base rules of flake8, the following flake8 extensions are also used:

- [flake8-simplify](https://github.com/MartinThoma/flake8-simplify) - for simplifying expressions.
- [flake8-comprehensions](https://github.com/adamchainz/flake8-comprehensions) - for list/set/dict comprehensions.
- [flake8-bugbear](https://github.com/PyCQA/flake8-bugbear) - for general bugs and design problems.
- [flake8-pyi](https://github.com/ambv/flake8-pyi) - for checking stub files.
- [flake8-docstrings](https://gitlab.com/pycqa/flake8-docstrings) - for checking docstrings. Currently not enforced
  fully.
- [flake8-pytest-style](https://github.com/m-burst/flake8-pytest-style) - for pytest.
- [flake8-eradicate](https://github.com/wemake-services/flake8-eradicate) - for removing commented out code.

## Continuous Integration

The CI that runs on every commit does the following:

- Perform static analysis using mypy and flake8. This is performed before tests, as it does not require dependencies to be installed.
- Unit testing using pytest. This checks the code is working as expected.
