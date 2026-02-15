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
python oscillator.py
```

This runs the full oscillator loop: Genesis4 undulator passes interleaved
with Fourier-optics cavity propagation. Runtime is approximately 60 s/pass
locally (~6.7 hours for 400 passes).

## Running on Koa HPC

The UH Koa cluster provides significantly faster execution with MPI
parallelism (~32 s/pass with 4 ranks).

### Single Pass with MPI

```bash
mpirun -np 4 genesis4 pass_0001.in
```

### Batch Job

Submit via SLURM. The oscillator script handles pass-to-pass file
management; each pass writes its output HDF5 file which the cavity
propagation code reads and modifies before the next pass.

## Output Files

Each pass produces an HDF5 output file (`pass_NNNN.out.h5`) containing:

- Radiation field on the transverse grid at the undulator exit
- Beam phase space
- Power along the bunch

The cavity propagation step modifies the radiation field and writes
it as the seed for the next Genesis4 pass.

## Performance Notes

| Configuration | Time per pass | Total (400 passes) |
|--------------|--------------|-------------------|
| Local (single core) | ~60 s | ~6.7 hours |
| Koa (4 MPI ranks) | ~32 s | ~3.6 hours |

Runtime is dominated by the Genesis4 undulator calculation. The
Fourier-optics cavity propagation is negligible in comparison.
