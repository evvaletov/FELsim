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

echo "C3v4d Job $SLURM_JOB_ID on $(hostname) at $(date)"
echo "COSY binary: $TMP_COSY"
echo "Python: $(python3 --version)"

# Approach C: Two-phase (stability-only → capped Twiss)
# v4 run reached instability=3.88 (approaching boundary at 1.0) but
# Phase 2 never activated because phase1_budget was too large.
# Fixed: phase1_budget = max_eval (one restart), then capped Twiss MSE.
echo ""
echo "========== Approach C: Two-phase objective =========="
python3 -u C3v4_cosy_mge_opt.py \
    --approach C \
    --fox-source input.fox \
    --sigma 2.0 \
    --max-eval 5000 \
    --restarts 15 \
    --save result_cosy_mge_v4C.json

# Approach D: Polish from best Approach C result
if [ -f result_cosy_mge_v4C.json ]; then
    echo ""
    echo "========== Approach D: Polish from Approach C =========="
    python3 -u C3v4_cosy_mge_opt.py \
        --approach D \
        --fox-source input.fox \
        --sigma 0.3 \
        --max-eval 5000 \
        --restarts 10 \
        --warm-start result_cosy_mge_v4C.json \
        --save result_cosy_mge_v4D.json
else
    echo ""
    echo "Skipping Approach D: result_cosy_mge_v4C.json not found"
fi

# Cleanup
rm -f "$TMP_COSY"
echo ""
echo "All approaches completed at $(date)"
