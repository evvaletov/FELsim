"""Generate chicane dipole fieldmap plots for documentation."""

import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
FIELDMAP_PATH = os.path.join(ROOT_DIR, 'fields', 'chicane_dipole_fieldmap.dat')
CSV_PATH = os.path.join(ROOT_DIR, 'fields', 'calculation', 'UH_chicane_fringe.csv')
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'images')

# Factor by which the old fieldmap was erroneously scaled
OLD_SCALE_FACTOR = 0.8351818473537908


def read_fieldmap(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[1].strip())
    deltas = float(lines[2].strip())
    by = np.array([float(lines[3 + i].strip()) for i in range(n)])
    z = (np.arange(n) - (n - 1) / 2) * deltas * 1000  # mm
    return z, by


def read_csv(path):
    z, by = [], []
    with open(path) as f:
        for row in csv.reader(f):
            z.append(float(row[0]))
            by.append(float(row[1]) / 10000)  # Gauss → T
    return np.array(z), np.array(by)


def compute_fwhm(z, by):
    half_max = by.max() / 2
    above = np.where(by >= half_max)[0]
    z_left = np.interp(half_max, by[above[0]-1:above[0]+1], z[above[0]-1:above[0]+1])
    z_right = np.interp(half_max, by[above[-1]:above[-1]+2][::-1], z[above[-1]:above[-1]+2][::-1])
    return z_right - z_left


def plot_profile(z_fm, by_fm, z_csv, by_csv):
    fig, ax = plt.subplots(figsize=(8, 4.5))

    ax.plot(z_fm, by_fm, 'b-', linewidth=1.2, label='Corrected fieldmap (201 pts)')
    ax.scatter(z_csv, by_csv, c='red', s=12, zorder=5, label='OPERA-3D source (132 pts)')

    peak = by_fm.max()
    fwhm = compute_fwhm(z_fm, by_fm)
    ax.annotate(f'Peak: {peak:.4f} T\nFWHM: {fwhm:.1f} mm',
                xy=(z_fm[np.argmax(by_fm)], peak), xytext=(30, peak * 0.75),
                fontsize=9, bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.8),
                arrowprops=dict(arrowstyle='->', color='gray'))

    ax.set_xlabel('z (mm)')
    ax.set_ylabel('$B_y$ (T)')
    ax.set_title('MkIII Chicane Dipole — Corrected Fieldmap Profile')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'fieldmap_profile.png'), dpi=150)
    plt.close(fig)
    print(f'  fieldmap_profile.png  (peak={peak:.4f} T, FWHM={fwhm:.1f} mm)')


def plot_correction(z_fm, by_fm):
    by_old = by_fm * OLD_SCALE_FACTOR

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(z_fm, by_old, 'r--', linewidth=1.0, label=f'Old (×{OLD_SCALE_FACTOR:.4f} scaling)')
    ax.plot(z_fm, by_fm, 'b-', linewidth=1.2, label='Corrected')

    # Shade the difference
    ax.fill_between(z_fm, by_old, by_fm, alpha=0.15, color='blue')

    peak_old = by_old.max()
    peak_new = by_fm.max()
    ax.annotate(f'Old peak: {peak_old:.4f} T\nNew peak: {peak_new:.4f} T\n'
                f'Δ: +{(peak_new - peak_old):.4f} T ({(peak_new/peak_old - 1)*100:.1f}%)',
                xy=(0, (peak_old + peak_new) / 2), xytext=(35, peak_new * 0.65),
                fontsize=9, bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', alpha=0.8),
                arrowprops=dict(arrowstyle='->', color='gray'))

    ax.set_xlabel('z (mm)')
    ax.set_ylabel('$B_y$ (T)')
    ax.set_title('MkIII Chicane Dipole — Before/After Fieldmap Correction')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'fieldmap_correction.png'), dpi=150)
    plt.close(fig)
    print(f'  fieldmap_correction.png  (old peak={peak_old:.4f} T, new peak={peak_new:.4f} T)')


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    z_fm, by_fm = read_fieldmap(FIELDMAP_PATH)
    z_csv, by_csv = read_csv(CSV_PATH)

    print('Generating fieldmap plots:')
    plot_profile(z_fm, by_fm, z_csv, by_csv)
    plot_correction(z_fm, by_fm)
    print('Done.')


if __name__ == '__main__':
    main()
