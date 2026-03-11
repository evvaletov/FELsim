#!/bin/bash
#SBATCH --job-name=cosy_mge_opt
#SBATCH --partition=shared-long
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=cosy_mge_opt_%j.log
#SBATCH --error=cosy_mge_opt_%j.err

cd /mnt/lustre/koa/scratch/evaletov/cosy_mge_opt

# Copy COSY binary to node-local /tmp (Lustre is noexec)
TMP_COSY="/tmp/cosy_${SLURM_JOB_ID}"
cp cosy "$TMP_COSY"
chmod +x "$TMP_COSY"
export COSY_BIN="$TMP_COSY"

# Skip broken user numpy, use system numpy + local cma
export PYTHONNOUSERSITE=1
export PYTHONPATH="$(pwd)"

echo "Job $SLURM_JOB_ID on $(hostname) at $(date)"
echo "COSY binary: $TMP_COSY"
echo "Python: $(python3 --version)"

python3 -u koa_cosy_mge_opt.py \
    --fox-source input.fox \
    --sigma 0.5 \
    --max-eval 8000 \
    --restarts 9 \
    --warm-start koa_cosy_mge_result.json \
    --save result_cosy_mge_v3.json

# Cleanup
rm -f "$TMP_COSY"
echo "Job completed at $(date)"
