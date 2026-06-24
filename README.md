# Chromatin_MSM
![workflow](image/workflow.png)

Python package combining two tools for chromatin fibres sampling:

- **`Chromatin_MSM.adaptive`** — adaptive MSM sampling pipeline
- **`Chromatin_MSM.feature_writer`** — PLUMED collective-variable feature generator

The MSM construction, validation and analysis can be found in `MSM_analysis/MSM_FMRC.ipynb`

## Install

```bash
pip install -e .
```

`OpenCGChromatin` is not on PyPI and is only required by `Chromatin_MSM.adaptive`
(the simulation/MSM pipeline, not the feature writer). Install it separately from
https://github.com/CollepardoLab/OpenCGChromatin and make sure it's importable
before running the adaptive-sampling. A working PLUMED installation (with the
Python bindings) and a CUDA-enabled OpenMM build are also required for `adaptive`.

### Required packages

- `OpenCGChromatin` ([GitHub, not PyPI/conda](https://github.com/CollepardoLab/OpenCGChromatin))
- `FMRC` ([GitHub, not PyPI/conda](https://github.com/Mingyuan00/Flow_Matching_for_RC))
- `python` (3.9 at least)
- `openmm` (8.1.2)
- `plumed` (2.9.2 / `py-plumed` 2.9.2)
- `mdtraj`
- `deeptime`
- `numpy`
- `pandas`
- `pytorch` (2.7.1 with CUDA 12.6 support)
- `matplotlib` 
- `scipy`
- `networkx`

## Usage

```bash
# Generate PLUMED CV definitions from a reference PDB
nucleosome-feature-writer --pdb ref/15/15rep.pdb --linker-length 15

# Run the adaptive MSM sampling pipeline
nucleosome-adaptive-sampling --n-rounds 8 --n-seed-per-round 16 --linker-length 15
```

Both are also usable as library functions, e.g. `from Chromatin_MSM.feature_writer import write_features`.

See `scripts/new.sh` for an example SLURM submission script for the adaptive-sampling
pipeline on a GPU cluster.

## Detailed docs

Full original documentation for each tool (directory layout, all CLI options,
feature descriptions, examples) is preserved under `docs/`:

- [`docs/adaptive.md`](docs/adaptive.md) — adaptive MSM sampling pipeline
- [`docs/feature_writer.md`](docs/feature_writer.md) — PLUMED feature generator

## Layout

```
Chromatin_MSM/
├── pyproject.toml
├── docs/
│   ├── adaptive.md              # adaptive sampling docs
│   ├── feature_writer.md        # feature generation docs
│   └── reference_structures.md  # ref/ docs
├── scripts/
│   └── new.sh                   # SLURM submission script
├── src/Chromatin_MSM/
│   ├── __init__.py
│   ├── adaptive.py               # adaptive MSM sampling pipeline
│   └── feature_writer.py         # PLUMED feature generator
├── ref/                          # reference structures, by linker length (bp)
│   ├── 15/15rep.pdb, core1-4.pdb, DNA.pdb, htail*_A-D.pdb, ...
│   ├── 20/ ... 22/ ... 25/ ... 25_A/ ... 25_AT/ ... 26/ ... 27/ ...
│   └── 28/ ... 29/ ... 30/ ... 30_A/ ... 30_AT/ ... 58/ ...
└── MSM_analysis/                 # MSM construction, validation and analysis notebooks
    └── MSM_FMRC.ipynb
```

Four functional areas:

- **Adaptive sampling** — `src/Chromatin_MSM/adaptive.py` ([docs](docs/adaptive.md))
- **Feature generation** — `src/Chromatin_MSM/feature_writer.py` ([docs](docs/feature_writer.md))
- **Reference structures** — `ref/` ([docs](docs/reference_structures.md))
- **MSM analysis** — `MSM_analysis/MSM_FMRC.ipynb`

Each module's original detailed documentation is preserved in its module
docstring/CLI `--help`.
