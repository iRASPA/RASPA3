"""Combine the ln Pi(N) of the TMMC windows in the directories 0, 1, ... and reweight it to other pressures.

Each window directory holds 'tmmc/tmmc_statistics.txt' written by RASPA3 ('MonteCarloTransitionMatrix'). The
windows are stitched on the macrostates they share (matched by N, the offset taken as the mean difference over the
whole overlap), the combined distribution is normalized, and its gas peak, liquid peak and the saddle between them
are located. The distribution is reweighted to other pressures with the Peng-Robinson fugacity of the pure
component, ln Pi(N; f1) = ln Pi(N; f0) + N ln(f1 / f0), not with the ratio of the pressures.

Run from this directory: python combine-tmmc-data.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema

# the state point of the example (simulation.json) and the critical constants of methane (methane.json)
TEMPERATURE = 111.58  # [K]
REFERENCE_PRESSURE = 36000.0  # [Pa], the 'ExternalPressure' every window ran at
CRITICAL_TEMPERATURE = 190.564  # [K]
CRITICAL_PRESSURE = 4599200.0  # [Pa]
ACENTRIC_FACTOR = 0.01142

STATISTICS_FILE = os.path.join("tmmc", "tmmc_statistics.txt")
FONT_SIZE = 15
LINE_WIDTH = 2


def peng_robinson_fugacity(temperature, pressure, tc=CRITICAL_TEMPERATURE, pc=CRITICAL_PRESSURE,
                           omega=ACENTRIC_FACTOR):
    """Fugacity [Pa] of a pure Peng-Robinson fluid, from the root of lowest Gibbs energy."""
    gas_constant = 8.314462618
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    alpha = (1.0 + kappa * (1.0 - np.sqrt(temperature / tc))) ** 2
    a = 0.45724 * gas_constant ** 2 * tc ** 2 / pc * alpha
    b = 0.07780 * gas_constant * tc / pc
    big_a = a * pressure / (gas_constant * temperature) ** 2
    big_b = b * pressure / (gas_constant * temperature)
    roots = np.roots([1.0, -(1.0 - big_b), big_a - 3.0 * big_b ** 2 - 2.0 * big_b,
                      -(big_a * big_b - big_b ** 2 - big_b ** 3)])
    roots = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > big_b]

    def ln_phi(z):
        sqrt2 = np.sqrt(2.0)
        return (z - 1.0 - np.log(z - big_b) - big_a / (2.0 * sqrt2 * big_b) *
                np.log((z + (1.0 + sqrt2) * big_b) / (z + (1.0 - sqrt2) * big_b)))

    # the stable root has the lowest ln phi (the lowest molar Gibbs energy at fixed T and P)
    return pressure * np.exp(min(ln_phi(z) for z in roots))


def read_window(directory):
    """ln Pi(N) of one window: the integer states (lambda bin 0 of an (N, lambda) chain)."""
    data = pd.read_csv(os.path.join(directory, STATISTICS_FILE), comment="#", sep=r"\s+")
    if "k" in data.columns:
        data = data[data["k"] == 0]
    return data[["N", "lnpi"]].reset_index(drop=True)


def combine_windows(directories):
    """Stitch the windows on the macrostates they share; each window is shifted onto the ones before it."""
    combined = None
    for directory in directories:
        window = read_window(directory)
        if combined is None:
            combined = window
            continue
        overlap = pd.merge(combined, window, on="N", suffixes=("_combined", "_window"))
        if overlap.empty:
            raise ValueError(f"window '{directory}' shares no macrostate with the windows before it")
        shift = np.mean(overlap["lnpi_combined"] - overlap["lnpi_window"])
        window = window.assign(lnpi=window["lnpi"] + shift)
        combined = pd.concat([combined, window[window["N"] > combined["N"].max()]], ignore_index=True)
    return combined.sort_values("N").reset_index(drop=True)


def normalized(lnpi):
    shifted = lnpi - np.max(lnpi)
    return shifted - np.log(np.sum(np.exp(shifted)))


def peaks_and_saddle(n, lnpi):
    """The gas peak, the liquid peak and the deepest saddle between them; for a unimodal ln Pi the single peak
    is both, and there is no saddle (None)."""
    lnpi = np.asarray(lnpi)
    n = np.asarray(n)
    maxima = argrelextrema(lnpi, np.greater)[0]
    if maxima.size == 0:
        maxima = np.array([int(np.argmax(lnpi))])
    gas, liquid = maxima.min(), maxima.max()
    minima = argrelextrema(lnpi, np.less)[0]
    saddles = minima[(minima > gas) & (minima < liquid)]
    saddle = saddles[np.argmin(lnpi[saddles])] if saddles.size > 0 else None
    return {
        "gas-N": n[gas], "gas-lnpi": lnpi[gas],
        "liq-N": n[liquid], "liq-lnpi": lnpi[liquid],
        "Saddle-N": n[saddle] if saddle is not None else np.nan,
        "Saddle-lnpi": lnpi[saddle] if saddle is not None else np.nan,
    }


def plot(n, lnpi, file_name, extremes):
    values = [extremes["gas-lnpi"], extremes["liq-lnpi"]]
    if not np.isnan(extremes["Saddle-lnpi"]):
        values.append(extremes["Saddle-lnpi"])
    figure = plt.figure()
    axis = figure.add_subplot(111)
    axis.plot(n, lnpi, color="orange", linewidth=LINE_WIDTH)
    axis.set_xlabel("Loading [Number of Molecules]", fontsize=FONT_SIZE, fontweight="bold")
    axis.set_ylabel("ln Pi", fontsize=FONT_SIZE, fontweight="bold")
    axis.set_ylim([min(values) - 10, max(values) + 10])
    figure.savefig(file_name, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main():
    directories = sorted((d for d in os.listdir(".") if d.isdigit() and os.path.exists(os.path.join(d, STATISTICS_FILE))),
                         key=int)
    if not directories:
        raise SystemExit(f"no window directories with '{STATISTICS_FILE}' found")

    data = combine_windows(directories)
    data["normalized-lnpi"] = normalized(data["lnpi"].to_numpy())

    reference_fugacity = peng_robinson_fugacity(TEMPERATURE, REFERENCE_PRESSURE)
    rows = []
    extremes = peaks_and_saddle(data["N"], data["normalized-lnpi"])
    rows.append({"Pressure": REFERENCE_PRESSURE, "Fugacity": reference_fugacity, **extremes})
    plot(data["N"], data["normalized-lnpi"], "combined_lnpi.png", extremes)

    for pressure in np.linspace(42500.0, 50000.0, 16):
        fugacity = peng_robinson_fugacity(TEMPERATURE, pressure)
        column = f"lnpi_{pressure:g}"
        data[column] = normalized(
            (data["normalized-lnpi"] + data["N"] * np.log(fugacity / reference_fugacity)).to_numpy())
        extremes = peaks_and_saddle(data["N"], data[column])
        rows.append({"Pressure": pressure, "Fugacity": fugacity, **extremes})
        plot(data["N"], data[column], f"lnpi_reweightedTo_{pressure:g}.png", extremes)

    data.to_csv("combined_lnpi.csv", index=False)
    pd.DataFrame(rows).to_csv("combined_min_max-lnpi.csv", index=False)


if __name__ == "__main__":
    main()
