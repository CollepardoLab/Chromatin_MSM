# Adaptive Sampling of MSMs for Tetra-nucleosome Arrays

Adaptive sampling of MSMs pipeline for chromatin simulations. Each round uses a Markov State Model (MSM) built from previous trajectories to seed new simulations in under-sampled regions of conformational space.


## Requirements
- python (> 3.9), OpenMM, OpenCGChromatin (https://github.com/CollepardoLab/OpenCGChromatin), DeepTime, MDTraj, PLUMED


## Directory Layout

```
.
├── adaptive.py          # main pipeline script
├── new.sh               # SLURM submission script
├── mdp_features.dat     # PLUMED input file defining CVs
├── ref/
│   └── 15rep.pdb        # reference topology
└── traj_dat/
    └── *.xtc            # original trajectories
```

Output directories:

```
CV/                  # PLUMED COLVAR files (all rounds)
adaptive_new1/       # round 1 simulation output
│   seed_0/
│   │   ├── trajectory.xtc
│   │   ├── seed_frame.pdb
│   │   └── state_data.out
│   seed_1/ ...
adaptive_new2/ ...
round1_seed.txt      # seed frame indices saved each round
```

---

## Running on a SLURM Cluster

Submit with:

```bash
sbatch new.sh
```

To pass custom parameters, edit the `python adaptive.py` line in `new.sh`:

```bash
python adaptive.py --n-rounds 5 --n-seed-per-round 10 --linker-length 15
```

---

## Command-Line Options

All parameters have sensible defaults and are optional.

### MSM / TICA parameters

| Parameters | Default | Description |
|------|---------|-------------|
| `--tica-lagtime` | `100` | TICA lag time |
| `--dim` | `None` | Number of TICA dimensions to keep (`None` = auto via `--var-cutoff`) |
| `--var-cutoff` | `0.95` | Cumulative variance cutoff for automatic TICA dimensionality selection |
| `--msm-lagtime` | `250` | MSM lag time |
| `--n-metastable-sets` | `20` | Number of macrostates for PCCA+ |

### Adaptive sampling control

| Paramters | Description |
|------|---------|-------------|
| `--n-rounds` | Total number of adaptive rounds to run |
| `--start-round` |  Round to start from (use `> 1` to resume a previous run) |
| `--n-seed-per-round` | New simulations launched per round |
| `--sim-idx` | Internal simulation counter offset |

### Example

```bash
python adaptive.py \
  --tica-lagtime 100 \
  --var-cutoff 0.95 \
  --msm-lagtime 250 \
  --n-metastable-sets 15 \
  --n-rounds 8 \
  --n-seed-per-round 16 \
  --linker-length 15 \
  --n-nucleosomes 4
```

### Resuming a run

If a job is interrupted after round 3, restart from round 4 without re-running earlier rounds:

```bash
python adaptive.py --start-round 4 --n-rounds 10
```

The pipeline will reload the existing COLVAR files and trajectories from rounds 1–3 automatically.

---

