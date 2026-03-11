#!/bin/bash
#SBATCH --job-name=c3v4_mge
#SBATCH --partition=shared-long
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=c3v4_mge_%j.log
#SBATCH --error=c3v4_mge_%j.err

cd /mnt/lustre/koa/scratch/evaletov/cosy_mge_opt

# Copy COSY binary to node-local /tmp (Lustre is noexec)
TMP_COSY="/tmp/cosy_${SLURM_JOB_ID}"
cp cosy "$TMP_COSY"
chmod +x "$TMP_COSY"
export COSY_BIN="$TMP_COSY"

export PYTHONNOUSERSITE=1
export PYTHONPATH="$(pwd)"

echo "C3v4 Job $SLURM_JOB_ID on $(hostname) at $(date)"
echo "COSY binary: $TMP_COSY"
echo "Python: $(python3 --version)"

# Approach A: Cold start + capped objective (most likely to succeed)
echo ""
echo "========== Approach A: Cold start + capped objective =========="
python3 -u C3v4_cosy_mge_opt.py \
    --approach A \
    --fox-source input.fox \
    --sigma 2.0 \
    --max-eval 5000 \
    --restarts 15 \
    --save result_cosy_mge_v4A.json

# Approach B: Warm start from v1 + capped objective
echo ""
echo "========== Approach B: Warm start + capped objective =========="
python3 -u C3v4_cosy_mge_opt.py \
    --approach B \
    --fox-source input.fox \
    --sigma 1.0 \
    --max-eval 5000 \
    --restarts 15 \
    --warm-start koa_cosy_mge_result.json \
    --save result_cosy_mge_v4B.json

# Approach C: Two-phase (stability-only → Twiss)
echo ""
echo "========== Approach C: Two-phase objective =========="
python3 -u C3v4_cosy_mge_opt.py \
    --approach C \
    --fox-source input.fox \
    --sigma 2.0 \
    --max-eval 5000 \
    --restarts 15 \
    --save result_cosy_mge_v4C.json

# Cleanup
rm -f "$TMP_COSY"
echo ""
echo "All approaches completed at $(date)"
