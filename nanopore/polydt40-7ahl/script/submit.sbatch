#!/bin/bash
#SBATCH --partition=dept_gpu
#SBATCH --job-name=ordered_pore
#SBATCH --exclude=g[020-023]
#SBATCH --constraint=2080Ti|Xp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --no-requeue
set -euo pipefail
ulimit -c 0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -u /net/galaxy/home/koes/zhangwen/ahl+3M-20-dA5/four-pore-translocation/ordered_translocation.py
