"""Chromatin_MSM: adaptive MSM sampling and PLUMED feature generation for tetranucleosome CG simulations.

Submodules are not imported eagerly here because `adaptive` requires the
optional, non-PyPI `OpenCGChromatin` and CUDA-enabled OpenMM dependencies.
Import them explicitly: `from Chromatin_MSM import feature_writer` or
`from Chromatin_MSM import adaptive`.
"""

