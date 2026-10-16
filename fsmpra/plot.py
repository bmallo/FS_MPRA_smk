"""Static scientific overview of pilot measurements; no significance encoding."""

from pathlib import Path

from .io import new_output, read_json


def plot_pilot(pilot_path, output, min_reads=10):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if min_reads < 1:
        raise ValueError("min_reads must be positive")
    pilot = read_json(pilot_path)
    wt = next(v for v in pilot["variants"] if v["class"] == "WT")
    snvs = sorted((v for v in pilot["variants"] if v["class"] == "SNV"), key=lambda v: (v["oligo_pos1"], v["alt"]))
    lo = pilot["region"]["oligo_start1"]
    length = len(pilot["region"]["sequence"])
    positions = np.arange(lo, lo + length)
    heatmap = np.full((len(snvs), length), np.nan)
    accessibility = np.full((4, length), np.nan)
    coverage = np.full((4, length), np.nan)
    for i, v in enumerate(snvs):
        coverage["ACGT".index(v["alt"]), v["oligo_pos1"] - lo] = v["n_reads"]
        if v["n_reads"] < min_reads:
            continue
        heatmap[i] = np.asarray(v["tracks"]["tf_small"]["delta_vs_WT"], dtype=float)
        accessibility["ACGT".index(v["alt"]), v["oligo_pos1"] - lo] = v["summary"]["msp_common_site_delta"]
    colors = plt.get_cmap("RdBu_r").copy()
    colors.set_bad("#e5e7eb")

    def limit(values):
        finite = np.abs(values[np.isfinite(values)])
        return max(.01, float(np.percentile(finite, 98))) if finite.size else .01

    fig = plt.figure(figsize=(15, 13), constrained_layout=True)
    grid = fig.add_gridspec(5, 2, height_ratios=[1, 2.6, .8, .8, 1])
    ax = fig.add_subplot(grid[0, :])
    for track, label, color in [("msp", "MSP coverage (accessibility proxy)", "#237a57"),
                                ("nucleosome", "Nucleosome coverage", "#9163b5"),
                                ("tf_small", "Small TF-track footprints", "#c77c21")]:
        ax.plot(positions, wt["tracks"][track]["occupancy"], label=label, color=color, lw=1.3)
    ax.set(title=f"WT landscape · {wt['n_reads']:,} sampled molecules", ylabel="Called occupancy", ylim=(0, 1))
    ax.legend(loc="upper right", fontsize=8, ncol=3)

    ax = fig.add_subplot(grid[1, :])
    lim = limit(heatmap)
    image = ax.imshow(heatmap, aspect="auto", interpolation="nearest", cmap=colors, vmin=-lim, vmax=lim,
                      extent=(lo - .5, lo + length - .5, len(snvs) - .5, -.5))
    ticks = list(range(0, len(snvs), max(1, len(snvs) // 6)))
    ax.set_yticks(ticks)
    ax.set_yticklabels([snvs[i]["source_variant"] for i in ticks], fontsize=8)
    ax.set(title="Small-footprint change vs WT · rows are individual SNVs", ylabel="Mutation (oligo coordinates)")
    fig.colorbar(image, ax=ax, label="Δ called occupancy", fraction=.025, extend="both")

    ax = fig.add_subplot(grid[2, :])
    lim = limit(accessibility)
    image = ax.imshow(accessibility, aspect="auto", interpolation="nearest", cmap=colors, vmin=-lim, vmax=lim,
                      extent=(lo - .5, lo + length - .5, 3.5, -.5))
    ax.set_yticks(range(4))
    ax.set_yticklabels(list("ACGT"))
    radius = pilot["method"]["common_site_exclusion_radius_bp"]
    ax.set(title=f"Region-wide MSP coverage change · common sites exclude mutation ±{radius} bp", ylabel="Alternate base")
    fig.colorbar(image, ax=ax, label="Δ MSP fraction", fraction=.025, extend="both")

    ax = fig.add_subplot(grid[3, :])
    image = ax.imshow(coverage, aspect="auto", interpolation="nearest", cmap="viridis",
                      extent=(lo - .5, lo + length - .5, 3.5, -.5))
    ax.set_yticks(range(4))
    ax.set_yticklabels(list("ACGT"))
    ax.set(title="Sampled molecules per allele", ylabel="Alternate base", xlabel="Oligo position")
    fig.colorbar(image, ax=ax, label="Molecules", fraction=.025)

    edges = np.asarray(pilot["method"]["histogram_edges"])
    for slot, key, title, xlabel in [
        (0, "methylation_histogram", "WT global methylation: full range retained", "Global m6A/AT (source metric)"),
        (1, "msp_fraction_histogram", "WT molecule-level accessibility proxy", "Fraction of mutable element covered by MSPs")]:
        ax = fig.add_subplot(grid[4, slot])
        vals = np.asarray(wt[key], dtype=float)
        if vals.sum():
            vals /= vals.sum()
        ax.bar(edges[:-1], vals, width=np.diff(edges), align="edge", color="#3d6e95", edgecolor="white")
        ax.set(title=title, xlabel=xlabel, ylabel="Fraction of measured WT", xlim=(0, 1))
    fraction = pilot["sampling"]["fraction"]
    fig.suptitle(f"{pilot['sample']['id']} · descriptive {fraction:.0%} pilot sample\n"
                 f"{length}-bp mutable element · effects shown for ≥{min_reads} molecules · no significance tests\n"
                 "Gray = reference base, missing data, or low coverage; effect color limits use the 98th percentile",
                 fontsize=12)
    out = new_output(output)
    fig.savefig(out / "overview.png", dpi=150)
    fig.savefig(out / "overview.pdf")
    plt.close(fig)
    return {"png": str(out / "overview.png"), "pdf": str(out / "overview.pdf")}
