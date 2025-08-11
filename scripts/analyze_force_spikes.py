#!/usr/bin/env python3
import argparse
import numpy as np


def read_pullf_xvg(path: str):
    times = []
    forces = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    times.append(float(parts[0]))
                    forces.append(float(parts[1]))
                except ValueError:
                    continue
    return np.array(times), np.array(forces)


def spike_score(forces: np.ndarray, zscore_thresh: float = 3.0):
    if forces.size == 0:
        return 1e9, 0, 0.0
    mu = np.mean(forces)
    sd = np.std(forces) if np.std(forces) > 1e-9 else 1.0
    z = np.abs((forces - mu) / sd)
    spikes = (z > zscore_thresh).sum()
    # Combine metrics: more spikes and higher mean are worse
    score = spikes * 1000.0 + float(mu)
    return score, spikes, float(mu)


def main():
    parser = argparse.ArgumentParser(description="Analyze force spikes across SMD windows")
    parser.add_argument('--inputs', nargs='+', required=True, help='List of pullf.xvg files')
    parser.add_argument('--output', required=True, help='Output summary text file')
    args = parser.parse_args()

    rows = []
    for path in args.inputs:
        t, f = read_pullf_xvg(path)
        score, nspikes, meanf = spike_score(f)
        rows.append((path, score, nspikes, meanf))

    rows.sort(key=lambda r: r[1])

    with open(args.output, 'w') as out:
        out.write("Window\tScore\tSpikes\tMeanForce\n")
        for path, score, nspikes, meanf in rows:
            out.write(f"{path}\t{score:.3f}\t{nspikes}\t{meanf:.3f}\n")
        if rows:
            out.write(f"\nBest window: {rows[0][0]} (score {rows[0][1]:.3f})\n")

    print(f"Wrote {args.output}")


if __name__ == '__main__':
    main() 