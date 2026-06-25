from OpenCGChromatin import *
from openmm import app, unit, Platform
import openmm as mm
import plumed
from deeptime.decomposition import TICA
from deeptime.clustering import KMeans
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.markov import TransitionCountEstimator
from deeptime.plots import plot_implied_timescales
from deeptime.util.validation import implied_timescales
from deeptime.plots.chapman_kolmogorov import plot_ck_test
import numpy as np
from numpy.random import multinomial
import pandas as pd
import os
import glob
import sys
import mdtraj as md
import multiprocessing
from multiprocessing import Pool
from copy import deepcopy
import subprocess
import argparse

colvar_prefix      = 'CV'
n_original_sims    = 320
feature_file       = 'mdp_features.dat'
traj_dir           = 'traj_dat'
# adaptive sampling parameters
tica_lagtime       = 100
dim                = None
var_cutoff         = 0.95
msm_lagtime        = 250
n_metastable_sets  = 20
n_rounds           = 10
sim_idx            = 0
start_round        = 1
n_seed_per_round   = 20
# simulation parameters
sim_steps          = 50000000 # 500 ns simulation length
save_steps         = 100000   # 100 ps saving interval
linker_length      = 15
N_nucleosomes      = 4
T_val              = 298
topology_file      = None  # derived from linker_length unless overridden; see default_topology_file()

def default_topology_file(linker_length):
    return f"ref/{linker_length}/{linker_length}rep.pdb"

def generate_actg_sequence(n):
    pattern = 'ACTG'
    return (pattern * (n // len(pattern) + 1))[:n]

def generate_even_dyads_sequence(N_nucleosomes, linker_length):
    # This is the W601 DNA sequence from literature
    w601 = ('ATCGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCT'
            'GTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCGAT')
    # Start with a short ACTG sequence for the half-linker on the left end.
    DNA_sequence = generate_actg_sequence(linker_length // 2)
    for _ in range(N_nucleosomes - 1):
        DNA_sequence += w601 + generate_actg_sequence(linker_length)
    DNA_sequence += w601 + generate_actg_sequence(linker_length // 2)
    return DNA_sequence

def process_trajectories(traj_glob, feature_file, out_dir, start_index=1):
    os.makedirs(out_dir, exist_ok=True)
    traj_files = sorted(glob.glob(traj_glob))
    if not traj_files:
        print(f"[WARNING] no trajectories match {traj_glob}")
        return
    for offset, traj in enumerate(traj_files):
        idx = start_index + offset
        cmd = ["plumed","driver","--plumed", feature_file, "--mf_xtc", traj]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            print(f"plumed driver failed:\n{res.stderr.strip()}")
            continue
        if not os.path.exists("COLVAR"):
            print(f"no COLVAR generated for {traj}")
            continue
        dst = os.path.join(out_dir, f"COLVAR_{idx}")
        os.replace("COLVAR", dst)

def read_features_indices(prefix_dir, start_idx, end_idx):
    data = []
    for i in range(start_idx, end_idx + 1):
        path = os.path.join(prefix_dir, f"COLVAR_{i}")
        try:
            df = plumed.read_as_pandas(path)
        except FileNotFoundError:
            print(f"[WARN] {path} not found, skipping")
            continue

        if "time" in df.columns:
            df = df.drop(columns=["time"])
        for col in list(df.columns):
            if col.startswith(("bend", "psi", "angle")):
            # drop raw angles columns
                df = df.drop(columns=[col])
        data.append(df.to_numpy(dtype=np.float32))

    return data

def run_TICA(data, data_supp, lagtime, dim=dim, var_cutoff=var_cutoff):
    data_all = data + data_supp
    tica = TICA(lagtime=lagtime, dim=dim, var_cutoff=var_cutoff)
    tica = tica.fit(data_all).fetch_model()
    tica_output = tica.transform(data)
    tica_output_concat = np.concatenate(tica_output)

    tica_output_supp = []
    for data_supp_i in data_supp:
        tica_output_supp_i = tica.transform(data_supp_i)
        tica_output_supp.append(tica_output_supp_i)

    # Include adaptive sampling data in the concatenated output
    if tica_output_supp:
        tica_output_supp_concat = np.concatenate(tica_output_supp)
        tica_output_concat = np.concatenate([tica_output_concat, tica_output_supp_concat])

    return tica, tica_output, tica_output_concat, tica_output_supp

def calculate_nmicro(tica_output_concat):
    # Heuristic approach to determine cluster number from htmd
    # https://github.com/Acellera/htmd/blob/master/htmd/adaptive/adaptivebandit.py
    n_microstates = int(max(100, np.round(0.6 * np.log10(tica_output_concat.shape[0] / 1000) * 1000 + 50)))
    return n_microstates

def run_kmeans(tica_output, tica_output_supp, tica_output_concat, n_microstates):
    clusters = KMeans(n_clusters=n_microstates, max_iter=500).fit_fetch(tica_output_concat)
    cluster_centers = clusters.cluster_centers
    # 1D array of cluster labels for each frame in the concatenated TICA output
    assignments_concat = clusters.transform(tica_output_concat)
    # 2D list of cluster labels contained time-series of each original trajectory-necessary for MSM construction
    assignments = assignments_concat.reshape(-1, tica_output.shape[1])

    assignments_supp = []
    for tica_output_supp_i in tica_output_supp:
        assignments_supp_i = clusters.transform(tica_output_supp_i)
        assignments_supp.append(assignments_supp_i)

    return assignments, assignments_concat, cluster_centers, assignments_supp

def build_MSM(msm_lagtime, assignments, assignments_supp):
    assignments_all = list(assignments) + assignments_supp
    counts = TransitionCountEstimator(lagtime=msm_lagtime, count_mode='sliding')
    counts = counts.fit_fetch(assignments_all)
    msm = MaximumLikelihoodMSM().fit_fetch(counts)
    return counts, msm

def check_connectivity(counts, n_microstates, connectivity_threshold=0., directed=True, min_largest_state_fraction=.95, max_n_components=1):
    sets = counts.connected_sets(connectivity_threshold=connectivity_threshold, directed=directed, sort_by_population=True)
    largest_component_fraction = len(sets[0]) / n_microstates
    ok = (largest_component_fraction >= min_largest_state_fraction) and (len(sets) <= max_n_components)
    return ok

def run_PCCA(msm, n_metastable_sets):
    pcca=msm.pcca(n_metastable_sets=n_metastable_sets)
    return pcca

def find_disconnected(counts, n_microstates, msm, n_metastable_sets, pcca):
    sets = counts.connected_sets(connectivity_threshold=0., directed=True, sort_by_population=True)
    disconnected = sets[1:]
    n_macro_disconnected = len(disconnected)
    disconnected_sets = {}
    # assign a unique macrostate label to each disconnected set
    for i, ds in enumerate(disconnected):
        macro_label = n_metastable_sets + i
        for state in ds:
            disconnected_sets[state] = macro_label
    pcca_assignments = np.zeros(n_microstates, dtype=int)
    stationary_distribution = np.zeros(n_microstates)
    connected_count = 0
    for micro in range(n_microstates):
        if micro in disconnected_sets:
            pcca_assignments[micro] = disconnected_sets[micro]
        else:
            pcca_assignments[micro] = pcca.assignments[connected_count]
            stationary_distribution[micro] = msm.stationary_distribution[connected_count]
            connected_count += 1
    return n_macro_disconnected, pcca_assignments, stationary_distribution

def count_macro(n_seed, pcca_assignments, assignments, assignments_concat, round_idx=None, assignments_supp=None):
    if round_idx is not None:
        print(f"[count_macro] Round {round_idx}")
    if assignments_supp:
        flat_supp = np.concatenate(assignments_supp, axis=0)
        assignments_concat = np.concatenate([assignments_concat, flat_supp], axis=0)
    # Two-level sampling: first sample macrostates, then microstates within each macrostate
    # macrostate-level
    # map macrostate assignment to simulation frames
    macro_timeseries = pcca_assignments[assignments_concat]
    # macrostate sort and count
    unique_macro, counts_macro = np.unique(macro_timeseries, return_counts=True)
    # set sampling probabilities inversely proportional to occupancy
    prob_macro = (1 / counts_macro) / np.sum(1 / counts_macro)
    # allocate no. seeds across macrostates according to prob_macro
    macrostate_seed = multinomial(n_seed, prob_macro)

    # microstate-level
    # microstate sort and count
    unique_micro, counts_micro = np.unique(assignments_concat, return_counts=True)
    seed_idx = []
    counts_micro_i_log = {}
    for macro_i, n_sample in enumerate(macrostate_seed):
        # locate the index of microstates not assigned to current selected macrostates
        not_macro_idx = np.where(pcca_assignments != unique_macro[macro_i])
        # let all entries corresponding to not_macro_idx = 0, therefore ignore them during selection
        counts_micro_i = deepcopy(counts_micro)
        counts_micro_i[not_macro_idx] = 0
        # let 1/0 = 0 (avoid division by zero warning)
        inverse_counts = np.divide(1, counts_micro_i, out=np.zeros_like(counts_micro_i, dtype=float), where=(counts_micro_i!=0))
        prob_micro_i = inverse_counts / np.sum(inverse_counts)
        microstate_seed = multinomial(n_sample, prob_micro_i)
        # Record selection statistics for visualization
        if n_sample != 0:
            macro_idx_log = unique_macro[macro_i]
            counts_micro_i_log[macro_idx_log] = [counts_micro_i,microstate_seed]
        for micro_i, n_sample in enumerate(microstate_seed):
            seed_idx = seed_idx + n_sample * [unique_micro[micro_i]]

    conf_seed = []
    for i, seed in enumerate(seed_idx):
        conf_idx = np.array(np.where(assignments==seed)).T
        conf_seed_frame = conf_idx[np.random.randint(conf_idx.shape[0], size=1), :][0]
        conf_seed.append(conf_seed_frame)

    return conf_seed

def initialize_cuda():
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
        return platform, {'Precision': 'Single'}
    except Exception as e:
        print("Error initializing CUDA:", e)
        return None, None

def run_one_adaptive_simulation(task):
    # run_params is passed explicitly (not read from module globals) because
    # multiprocessing's 'spawn' context re-imports this module fresh in each
    # worker, so any values main() set at runtime would otherwise be invisible
    # here and silently fall back to the module-level defaults.
    (traj_path, frame), round_idx, seed_index, gpu_id, run_params = task
    topology_file = run_params["topology_file"]
    linker_length = run_params["linker_length"]
    N_nucleosomes = run_params["N_nucleosomes"]
    sim_steps = run_params["sim_steps"]
    save_steps = run_params["save_steps"]
    T_val = run_params["T_val"]
    try:
        platform, properties = initialize_cuda()
        if platform is None:
            raise RuntimeError("CUDA platform not available")
        properties['CudaDeviceIndex'] = str(gpu_id)
        # Load directly from the path provided in the task (no globals / indexing)
        traj = md.load(traj_path, top=topology_file)
        seed_dir = os.path.join(f"adaptive_new{round_idx}", f"seed_{seed_index}")
        os.makedirs(seed_dir, exist_ok=True)
        traj[frame].save_pdb(os.path.join(seed_dir, "seed_frame.pdb"))
        pdb = app.PDBFile(os.path.join(seed_dir, "seed_frame.pdb"))
        positions = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        DNA_sequence = generate_even_dyads_sequence(N_nucleosomes, linker_length)
        nucleosome_sequence = N_nucleosomes * ['1kx5']
        first_dyad = linker_length//2 + 73
        dyad_positions = [first_dyad + i*(147+linker_length) for i in range(N_nucleosomes)]
        nucleosome_array = NucleosomeArray(nucleosome_sequence, DNA_sequence, dyad_positions)
        topology= nucleosome_array.topology
        system = get_system(positions, topology, {nucleosome_array.chain_id: nucleosome_array.globular_indices}, [nucleosome_array.dyad_positions], debye_length=0.8, constraints='breathing', periodic=False)
        integrator = mm.LangevinMiddleIntegrator(T_val*unit.kelvin, 0.01/unit.picosecond, 10*unit.femtosecond)
        simulation = app.Simulation(topology, system, integrator, platform=platform, platformProperties=properties)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(app.XTCReporter(os.path.join(seed_dir, "trajectory.xtc"), save_steps))
        simulation.reporters.append(app.StateDataReporter(os.path.join(seed_dir, "state_data.out"), save_steps, step=True, potentialEnergy=True, temperature=True, separator='\t'))
        simulation.step(sim_steps)
        return True
    except Exception as e:
        print(f"Error in round {round_idx}, seed {seed_index}: {e}")
        return False

def run_adaptive_parallel(tasks, max_concurrent=n_seed_per_round):
    ctx = multiprocessing.get_context('spawn')
    with ctx.Pool(processes=max_concurrent) as pool:
        return pool.map(run_one_adaptive_simulation, tasks)

def preflight_check(traj_dir, topology_file, feature_file, start_round):
    """Validate required inputs exist before doing any TICA/MSM work.

    Fails loudly with an actionable message instead of running for a while
    and dying deep inside deeptime (e.g. AssertionError: Data is empty.).
    """
    problems = []

    traj_files = glob.glob(os.path.join(traj_dir, "*.xtc"))
    if not traj_files:
        problems.append(
            f"No trajectories found in '{traj_dir}/*.xtc'. The original seed "
            f"trajectories must be placed there before running."
        )

    if not os.path.isfile(topology_file):
        problems.append(
            f"Topology file '{topology_file}' not found. Either place a PDB "
            f"there, or pass --topology-file pointing at the right structure "
            f"for --linker-length."
        )

    if not os.path.isfile(feature_file):
        problems.append(
            f"Feature file '{feature_file}' not found. Generate it first with: "
            f"nucleosome-feature-writer --pdb {topology_file} "
            f"--linker-length <N> --output {feature_file}"
        )

    if problems:
        print("\nPreflight check failed:\n", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        print(
            "\nSee scripts/setup_adaptive_run.sh for a helper that stages "
            "the topology symlink and feature file for a run directory.\n",
            file=sys.stderr,
        )
        sys.exit(1)

def main():
    global tica_lagtime, dim, var_cutoff, msm_lagtime, n_metastable_sets
    global n_rounds, sim_idx, start_round, n_seed_per_round
    global linker_length, N_nucleosomes, topology_file, feature_file, traj_dir

    parser = argparse.ArgumentParser(description="Adaptive sampling pipeline")
    parser.add_argument("--tica-lagtime",      type=int,   default=tica_lagtime)
    parser.add_argument("--dim",               type=int,   default=dim)
    parser.add_argument("--var-cutoff",        type=float, default=var_cutoff)
    parser.add_argument("--msm-lagtime",       type=int,   default=msm_lagtime)
    parser.add_argument("--n-metastable-sets", type=int,   default=n_metastable_sets)
    parser.add_argument("--n-rounds",          type=int,   default=n_rounds)
    parser.add_argument("--sim-idx",           type=int,   default=sim_idx)
    parser.add_argument("--start-round",       type=int,   default=start_round)
    parser.add_argument("--n-seed-per-round",  type=int,   default=n_seed_per_round)
    parser.add_argument("--linker-length",     type=int,   default=linker_length)
    parser.add_argument("--n-nucleosomes",     type=int,   default=N_nucleosomes)
    parser.add_argument("--traj-dir",          type=str,   default=traj_dir,
                         help="Directory containing the original seed trajectories (*.xtc)")
    parser.add_argument("--topology-file",     type=str,   default=None,
                         help="Reference topology PDB. Default: ref/<linker-length>/<linker-length>rep.pdb")
    parser.add_argument("--feature-file",      type=str,   default=feature_file,
                         help="PLUMED feature definition file (from nucleosome-feature-writer)")
    args = parser.parse_args()

    tica_lagtime      = args.tica_lagtime
    dim               = args.dim
    var_cutoff        = args.var_cutoff
    msm_lagtime       = args.msm_lagtime
    n_metastable_sets = args.n_metastable_sets
    n_rounds          = args.n_rounds
    sim_idx           = args.sim_idx
    start_round       = args.start_round
    n_seed_per_round  = args.n_seed_per_round
    linker_length     = args.linker_length
    N_nucleosomes     = args.n_nucleosomes
    traj_dir          = args.traj_dir
    feature_file      = args.feature_file
    topology_file     = args.topology_file or default_topology_file(linker_length)

    preflight_check(traj_dir, topology_file, feature_file, start_round)

    # 1) Process original trajectories
    print("Processing original trajectories...")
    process_trajectories(os.path.join(traj_dir, "*.xtc"), feature_file, colvar_prefix, start_index=1)
    # load
    orig_paths = sorted(glob.glob(os.path.join(traj_dir, "*.xtc")))
    orig_trajs = [md.load(p, top=topology_file) for p in orig_paths]
    new_paths = []

    # 2) Initial MSM
    print("\nBuilding initial MSM...")
    prev_adaptive = (start_round - 1) * n_seed_per_round
    if prev_adaptive > 0:
        first_prev=n_original_sims+1
        last_prev=n_original_sims+prev_adaptive
        data_prev=read_features_indices(colvar_prefix,first_prev,last_prev)
    else:
        data_prev=[]
    data = read_features_indices(colvar_prefix, 1, n_original_sims)
    data_supp = data_prev
    tica, tica_output, tica_output_concat, tica_output_supp = run_TICA(data, data_supp, tica_lagtime, dim, var_cutoff)
    print(tica_output_concat.shape[0])
    n_microstates = calculate_nmicro(tica_output_concat)
    print(f"Estimated number of microstates: {n_microstates}")
    assignments, assignments_concat, cluster_centers, assignments_supp = run_kmeans(tica_output, tica_output_supp, tica_output_concat, n_microstates)
    print(assignments_concat.shape)
    counts, msm = build_MSM(msm_lagtime, assignments, assignments_supp)
    print("MSM built with", msm.transition_matrix.shape[0], "states.")
    # 3) Check connectivity
    connectivity_ok = check_connectivity(counts, n_microstates, connectivity_threshold=0., directed=True, min_largest_state_fraction=.95, max_n_components=1)
    if connectivity_ok:
        print('MSM is well-connected, skipping adaptive sampling.')
        sys.exit(0)
    pcca = run_PCCA(msm, n_metastable_sets)
    n_macro_disconnected, pcca_assignments, stationary_distribution = find_disconnected(counts, n_microstates, msm, n_metastable_sets, pcca)
    conf_seed = count_macro(n_seed_per_round, pcca_assignments, assignments, assignments_concat, sim_idx, assignments_supp=assignments_supp)
    seed_positions = conf_seed
    print(f"Initial seeding complete: {len(seed_positions)} seeds generated.")

    # 4) Adaptive rounds
    all_data = data.copy()
    all_data.extend(data_prev)
    # gather paths from rounds 1..start_round-1
    new_paths = []
    if start_round > 1:
        for r in range(1, start_round):
            new_paths.extend(sorted(glob.glob(f"adaptive_new{r}/seed_*/trajectory.xtc")))
    sim_idx = start_round - 1
    for round_i in range(start_round, n_rounds+1):
        print(f"\n=== ROUND {round_i} ===")
        np.savetxt(f"round{round_i}_seed.txt", seed_positions, fmt="%d")
        sim_idx += 1
        combined_paths = orig_paths + new_paths
        print(f"[resume] combined_paths: {len(combined_paths)} = "
      f"{len(orig_paths)} + {len(new_paths)}")
        # a) run sims
        run_params = {
            "topology_file": topology_file,
            "linker_length": linker_length,
            "N_nucleosomes": N_nucleosomes,
            "sim_steps": sim_steps,
            "save_steps": save_steps,
            "T_val": T_val,
        }
        tasks = [((combined_paths[i], f), round_i, k, k % 4, run_params)
             for k, (i, f) in enumerate(seed_positions)]
        run_adaptive_parallel(tasks, max_concurrent=n_seed_per_round)
        # b) process new CVs
        start_idx = n_original_sims + (round_i-1)*n_seed_per_round + 1
        process_trajectories(f"adaptive_new{round_i}/seed_*/trajectory.xtc",feature_file, colvar_prefix, start_idx)
        # c) load new data
        new_paths_this = sorted(glob.glob(f"adaptive_new{round_i}/seed_*/trajectory.xtc"))
        new_paths.extend(new_paths_this)
        new_data = read_features_indices(colvar_prefix,start_idx,start_idx+n_seed_per_round-1)
        all_data.extend(new_data)
        # d) refit MSM
        data_supp = all_data[n_original_sims:]
        tica, tica_output, tica_output_concat, tica_output_supp = run_TICA(data, data_supp, tica_lagtime, dim, var_cutoff)
        print(tica_output_concat.shape[0])
        n_microstates = calculate_nmicro(tica_output_concat)
        print(f"Estimated number of microstates: {n_microstates}")
        assignments, assignments_concat, cluster_centers, assignments_supp = run_kmeans(tica_output, tica_output_supp, tica_output_concat, n_microstates)
        counts, msm = build_MSM(msm_lagtime, assignments, assignments_supp)
        print("MSM built with", msm.transition_matrix.shape[0], "states.")
        pcca = run_PCCA(msm, min(n_metastable_sets, msm.transition_matrix.shape[0]))
        n_macro_disconnected, pcca_assignments, _ = find_disconnected(counts, n_microstates, msm, n_metastable_sets, pcca)
        # e) next seeds
        if round_i <= n_rounds:
            seed_positions = count_macro(n_seed_per_round, pcca_assignments, assignments, assignments_concat, round_i, assignments_supp=assignments_supp)
    print("\nAdaptive sampling complete :)")


if __name__ == "__main__":
    main()
