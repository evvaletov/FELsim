# Running Genesis4 Simulations

## Prerequisites

- Genesis4 compiled and on `$PATH`
- Python with numpy, h5py, scipy, matplotlib
- For HPC runs: MPI (mpirun)

## Local Execution

### Single Steady-State Pass

```bash
cd /home/evaletov/UH/GENESIS/UHFEL_undulator/genesis4
genesis4 uhfel_ss.in
```

### Multi-Pass Oscillator

```bash
cd /home/evaletov/UH/GENESIS/UHFEL_undulator/genesis4

# Default (SAMPLE=3, normal saturation)
python oscillator.py 400

# Superradiant regime (SAMPLE=1) — recommended
python oscillator.py 400 --sample 1

# With desynchronization
python oscillator.py 400 --sample 1 --desync 0.5

# Seeded start (10 W)
python oscillator.py 300 --seed 10

# Full options
python oscillator.py 400 --sample 1 --slen 0.01 --np 4 --lambda0 3.29e-6
```

### CLI Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `npass` | 400 | Number of cavity passes |
| `--np` | 1 | MPI ranks for Genesis4 |
| `--sample` | 3 | Wavelengths per time slice (**use 1 for superradiance**) |
| `--slen` | 0.0025 | Time window length in meters |
| `--desync` | 0.0 | Desynchronization parameter $d$ |
| `--seed` | 0 | Seed power in watts (0 = shot noise) |
| `--lambda0` | 3.29e-6 | Carrier wavelength |
| `--fresnel-kernel` | off | Use Fresnel instead of exact Helmholtz |
| `--aperture` | off | Enable soft mirror aperture |

## Running on Koa HPC

The UH Koa cluster provides significantly faster execution with MPI
parallelism.

### SLURM Script Template

```bash
#!/bin/bash
#SBATCH --job-name=uhfel-osc
#SBATCH --partition=shared
#SBATCH --account=uh
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --output=uhfel_%j.out

module load lang/Anaconda3
eval "$(conda shell.bash hook)"
conda activate scientific
export PATH="$HOME/.conda/envs/genesis4/bin:$PATH"

# IMPORTANT: unset SLURM variables that conflict with mpirun
unset SLURM_CPUS_PER_TASK
unset SLURM_TRES_PER_TASK

python oscillator.py 400 --sample 1 --np $SLURM_NTASKS
```

### Important Notes

- Genesis4 binary: `~/.conda/envs/genesis4/bin/genesis4` (v4.6.11)
- Python environment: `scientific` (has h5py, numpy, scipy)
- **Must unset** `SLURM_CPUS_PER_TASK` and `SLURM_TRES_PER_TASK` before `mpirun`
- Account: `uh` (not `hpc_evaletov1`)
- Sandbox partition: ~10 tasks/node max, 4h limit (for quick tests)
- Shared partition: up to 36h, recommended for production runs

## Output Files

Each pass produces an HDF5 output file (`pass_NNNN.out.h5`) containing:

- Radiation field on the transverse grid at the undulator exit
- Beam phase space
- Power along the bunch

The cavity propagation step modifies the radiation field and writes
it as the seed for the next Genesis4 pass. Intermediate files are
cleaned up automatically; only the final cavity field and
`oscillator_results.npz` are retained.

### Results File

`oscillator_results.npz` contains:

| Array | Shape | Description |
|-------|-------|-------------|
| `peak_powers` | `(npass,)` | Peak power per pass [W] |
| `pulse_energies` | `(npass,)` | Integrated pulse energy [J] |
| `power_profiles` | `(npass, nslices)` | Full temporal power profile |
| `rms_sizes` | `(npass,)` | RMS transverse spot size at peak [m] |
| `peak_slice_idx` | `(npass,)` | Index of peak slice |

Scalars: `npass`, `desync_d`, `ds`, `slen`, `sample`, `lambda0`.

### Plotting

```bash
python plot_diagnostics.py
```

Generates comparison plots in `results/plots/`: power evolution, pulse
profiles, duration, waterfall, gain per pass, RMS size, energy evolution.

## Performance Notes

| Configuration | SAMPLE | Time/pass | Total (400 passes) |
|--------------|--------|-----------|-------------------|
| Local (1 core) | 3 | ~60 s | ~6.7 hours |
| Koa (4 MPI) | 3 | ~32 s | ~3.6 hours |
| Koa (4 MPI) | 1 | ~100 s | ~11 hours |
| Koa (8 MPI) | 1 + SLEN=10mm | ~265 s | ~30 hours |

`SAMPLE=1` is ~3× slower than `SAMPLE=3` due to 3× more time slices.
Larger `SLEN` increases memory and runtime proportionally. Runtime is
dominated by the Genesis4 undulator calculation; Fourier-optics cavity
propagation is negligible.
