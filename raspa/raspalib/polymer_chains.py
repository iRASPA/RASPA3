"""Combines the chain-length resolved polymer chain statistics of independent polymerization runs.

Every run of a polymerization with 'ComputePolymerChains' (or a ``PropertyPolymerChains`` given to a system of
the Python driver) holds per key (chain length, topology, living, conversion bin) the number of chain samples and
the sums and sums of squares of the per-chain values (Rg^2, Rg, the eigenvalues of the gyration tensor, R^2, ...).
The runs are combined by adding these sums. The samples of one run are correlated (the same chains, sampled again
and again), so the error of a mean is taken from the scatter of the means of the runs about the combined mean,
each run weighted by its number of samples:

    mean = sum_r S_r / sum_r n_r
    error^2 = sum_r (n_r / n_bar)^2 (m_r - mean)^2 / (R (R - 1)),     m_r = S_r / n_r,  n_bar = sum_r n_r / R

with R the number of runs that have samples of the key (the error is NaN when fewer than two have).

    from raspalib import polymer_chains
    rows = polymer_chains.combine(glob.glob("run*/polymer_chains/polymer_chains.s0.json"))
    for row in rows:
        print(row["length"], row["topology"], row["mean"]["Rg^2"], row["error"]["Rg^2"])
"""

import json
import math
import os
from collections.abc import Iterable, Sequence

KEY_FIELDS = ("length", "topology", "living", "conversionBin")


def _load(source):
    """The JSON content of a run: a path, an already parsed dictionary, or a PropertyPolymerChains."""
    if isinstance(source, (str, os.PathLike)):
        with open(source) as stream:
            return json.load(stream)
    if isinstance(source, dict):
        return source
    if hasattr(source, "entries"):
        return {
            "moments": list(source.moments),
            "settings": {"conversionBins": source.conversion_bins},
            "entries": source.entries(),
        }
    raise TypeError(f"polymer_chains: cannot read a run from {type(source).__name__}")


def combine(sources: Iterable, by: Sequence[str] = ("length", "topology")) -> list[dict]:
    """Combines runs per key ``by`` (a subset of 'length', 'topology', 'living', 'conversionBin').

    Returns one dictionary per key, sorted, with the key fields, 'runs' (the number of runs with samples),
    'count' (the samples), and per value name 'mean', 'error' (from the scatter between the runs; NaN when fewer
    than two runs have samples of the key) and 'standard_deviation' (the width of the distribution over all samples), plus
    'end_to_end_ratio' <R^2>/<Rg^2>.
    """
    unknown = [field for field in by if field not in KEY_FIELDS]
    if unknown:
        raise ValueError(f"polymer_chains: unknown key fields {unknown}; choose from {KEY_FIELDS}")

    runs = [_load(source) for source in sources]
    if not runs:
        return []
    moments = runs[0]["moments"]
    for run in runs[1:]:
        if run["moments"] != moments:
            raise ValueError("polymer_chains: the runs accumulate different values")
    if "conversionBin" in by:
        bins = {run.get("settings", {}).get("conversionBins") for run in runs}
        if len(bins) > 1:
            # (a run that does not state its bins is listed as None, which does not sort with numbers)
            raise ValueError(
                f"polymer_chains: the runs use different numbers of conversion bins {sorted(bins, key=repr)}"
            )

    # per key, per run: count, sums and sums of squares
    width = len(moments)
    per_key: dict[tuple, list[list]] = {}
    for r, run in enumerate(runs):
        for entry in run["entries"]:
            key = tuple(entry[field] for field in by)
            slots = per_key.setdefault(key, [[0, [0.0] * width, [0.0] * width] for _ in runs])
            slot = slots[r]
            slot[0] += entry["count"]
            for m in range(width):
                slot[1][m] += entry["sum"][m]
                slot[2][m] += entry["sumOfSquares"][m]

    rows = []
    for key in sorted(per_key):
        slots = per_key[key]
        count = sum(slot[0] for slot in slots)
        if count == 0:
            continue
        sampled = [slot for slot in slots if slot[0] > 0]
        number_of_runs = len(sampled)
        mean_count = count / number_of_runs
        row = dict(zip(by, key))
        row["runs"] = number_of_runs
        row["count"] = count
        row["mean"], row["error"], row["standard_deviation"] = {}, {}, {}
        for m, name in enumerate(moments):
            total = sum(slot[1][m] for slot in slots)
            squares = sum(slot[2][m] for slot in slots)
            mean = total / count
            if number_of_runs > 1:
                # The mean of run r, of n_r samples, scatters as sigma^2 / n_r about the common mean: with the
                # weights w_r = n_r / N its error is sum_r w_r (m_r - mean)^2 / (R - 1), the classic error of a mean
                # for runs of equal length. (The squared weights over R (R - 1) took the scatter of the long runs as
                # that of all, and went to zero for runs of very different length.)
                scatter = sum((slot[0] / count) * (slot[1][m] / slot[0] - mean) ** 2 for slot in sampled)
                error = math.sqrt(scatter / (number_of_runs - 1))
            else:
                error = math.nan
            variance = (squares - count * mean * mean) / (count - 1) if count > 1 else 0.0
            row["mean"][name] = mean
            row["error"][name] = error
            row["standard_deviation"][name] = math.sqrt(max(variance, 0.0))
        rg2 = row["mean"].get("Rg^2", 0.0)
        row["end_to_end_ratio"] = row["mean"].get("R^2", 0.0) / rg2 if rg2 > 0.0 else math.nan
        rows.append(row)
    return rows


def write(rows: Sequence[dict], path: str, names: Sequence[str] = ("Rg^2", "Rg", "R^2", "k^2")) -> None:
    """Writes combined rows as a column table: the key fields, runs, samples, and per name the mean and error."""
    if not rows:
        with open(path, "w") as stream:
            stream.write("# no samples\n")
        return
    key_fields = [field for field in KEY_FIELDS if field in rows[0]]
    with open(path, "w") as stream:
        columns = key_fields + ["runs", "count"]
        for name in names:
            columns += [f"<{name}>", f"error({name})"]
        columns.append("<R^2>/<Rg^2>")
        for i, column in enumerate(columns, start=1):
            stream.write(f"# column {i}: {column}\n")
        for row in rows:
            fields = [str(row[field]) for field in key_fields] + [str(row["runs"]), str(row["count"])]
            for name in names:
                fields += [f"{row['mean'][name]:.6g}", f"{row['error'][name]:.6g}"]
            fields.append(f"{row['end_to_end_ratio']:.6g}")
            stream.write(" ".join(fields) + "\n")
