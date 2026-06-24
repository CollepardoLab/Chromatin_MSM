#!/bin/bash
#SBATCH --job-name=15_adp
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:4
#SBATCH --time=48:00:00
#SBATCH --account=ehpc330
#SBATCH --qos=acc_ehpc
#SBATCH -D .
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err

export SLURM_CPU_BIND=none
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

# Env
source ~/.bashrc
conda activate openmmtools
# --- MPS setup (per node, per job) ---
MPS_ROOT="$TMPDIR/mps-$USER-${SLURM_JOB_ID:-$$}-$(hostname)"
export CUDA_MPS_PIPE_DIRECTORY="$MPS_ROOT/pipe"
export CUDA_MPS_LOG_DIRECTORY="$MPS_ROOT/log"

# 3) Create and lock down permissions (don’t use $HOME on NFS)
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
chmod 700 "$MPS_ROOT" "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"

# Start MPS control daemon
nvidia-cuda-mps-control -d

# If you're just running one program, this is fine:
nucleosome-adaptive-sampling

# If you want multiple clients per GPU using MPS, e.g. two procs per GPU:
# for i in {1..8}; do
#   CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=50 srun --exclusive -N1 -n1 nucleosome-adaptive-sampling &
# done
# wait

# Stop MPS and clean up
echo quit | nvidia-cuda-mps-control
rm -rf "$MPS_ROOT"

