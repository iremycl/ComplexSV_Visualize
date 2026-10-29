"""Plot structural variants and copy number for one or more chromosomes.

The main plotting helpers are kept here as ordinary Python functions so they
can be reused from another script or imported in a notebook. The command-line
entry point at the bottom reads three tab-separated files and writes a
multi-page PDF.

Input files:

* coverage: chromosome, start, end, total copy number
* segmentation: chromosome, start, end, total copy number, minimum copy number
* BEDPE: structural-variant records with cluster and type columns in the
  positions expected by ``load_bedpe``
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from collections import defaultdict
import matplotlib.transforms as mtransforms

# ====== constants ======
MISSING = {"", ".", "-", "NA", None}
SV_MARKER = {"td": "^", "del": "v", "inv": "s", "inter": "D", "complex": "o"}
SV_MARKER_SIZE = 7
CN_YLIM = (0, 12)
CLUSTER_ALPHA_SV = 0.18    # shaded span alpha on SV panel
CLUSTER_ALPHA_CN = 0.10    # shaded span alpha on CN panel
INTER_COLOR = "deeppink"   # inter-chrom vertical guide + label color

# ====== Helper functions ======
def sort_chromosomes(chr_list):
    """Sort chromosomes like chr1..chr22, chrX, chrY, then others."""
    def key(c):
        c = str(c).replace("chr", "")
        if c.isdigit(): return (0, int(c))
        if c in {"X", "Y"}: return (1, {"X": 23, "Y": 24}[c])
        return (2, c)
    return sorted(chr_list, key=key)

def _sv_category(ctype: str) -> str:
    """Bucket SV type from column 18 into td/del/inv/inter/complex."""
    c = (ctype or "").lower()
    if c.startswith("td") or c.startswith("dup"): return "td"
    if c.startswith("del"):                       return "del"
    if c.startswith("inv"):                       return "inv"
    if c.startswith("inter"):                     return "inter"
    return "complex"  # everything else

# ====== parsing helpers ======
def load_coverage(path: str) -> pd.DataFrame:
    """Read binned coverage: chr, start, end, total_cn."""
    return pd.read_csv(
        path, sep="\t", header=None,
        names=["chr", "start", "end", "total_cn"],
        dtype={"chr": str, "start": np.int32, "end": np.int32, "total_cn": np.float64}
    )

def load_segmentation(path: str) -> pd.DataFrame:
    """Read segmentation: chr, start, end, total_cn, min_cn."""
    return pd.read_csv(
        path, sep="\t", header=None,
        names=["chr", "start", "end", "total_cn", "min_cn"],
        dtype={"chr": str, "start": np.int32, "end": np.int32,
               "total_cn": np.float64, "min_cn": np.float64}
    )

def load_bedpe(path: str) -> pd.DataFrame:
    """
    Read BEDPE + light normalization:
      - rename 0-based columns to friendlier names
      - derive cluster_short (suffix after last ':')
      - compute cluster_size and sv_type (simple/clustered/complex)
    """
    # IMPORTANT: indices here are 0-based (header=None).
    bedpe_names = {
        "0": "chr1", 
        "1": "start1", 
        "2": "end1",
        "3": "chr2", 
        "4": "start2", 
        "5": "end2",
        "6": "id", 
        "7": "qual", 
        "8": "strand1", 
        "9": "strand2",
        "12": "cluster_id",      # 13th col (1-based) in bedpe = 12 here
        "17": "complex_type",    # 18th col
        "18": "footprint_low",   # 19th col (optional)
        "19": "footprint_high",  # 20th col (optional)
    }

    bedpe_dtype = {
        0: str, 
        1: np.int64, 
        2: np.int64,
        3: str, 
        4: np.int64, 
        5: np.int64,
        6: str, 
        7: np.float64, 
        8: str, 
        9: str,
        10: str, 
        12: str, 
        17: str, 
        18: str, 
        19: str,
    }

    bedpe = (
        pd.read_csv(path, sep="\t", header=None, dtype=bedpe_dtype, low_memory=False)
          .rename(columns=str)
          .rename(columns=bedpe_names)
          .assign(
              # Short cluster label (after last ':'), keep as-is for NA-like tokens
              cluster_short=lambda df: df["cluster_id"].where(~df["cluster_id"].isin(MISSING))
                                         .str.rsplit(":", n=1).str[-1]
          )
          .sort_values(["chr1", "start1", "chr2", "start2"]).reset_index(drop=True)
    )

    # cluster_size across the whole file by short id (0 if NA)
    counts = bedpe["cluster_short"].value_counts()
    bedpe["cluster_size"] = bedpe["cluster_short"].map(counts).fillna(0).astype(int)

    # sv_type init
    bedpe["sv_type"] = "simple"
    # complex if complex_type starts with 'complex'
    bedpe.loc[bedpe["complex_type"].str.startswith("complex", na=False), "sv_type"] = "complex"
    # clustered if not complex and cluster_size >= 2
    bedpe.loc[(bedpe["sv_type"] != "complex") & (bedpe["cluster_size"] >= 2), "sv_type"] = "clustered"

    # quick sanity
    print("SV types:", bedpe["sv_type"].value_counts().to_dict())
    return bedpe

def load_data(coverage_file, segmentation_file, bedpe_file):
    """Convenience wrapper to read all three inputs."""
    return (
        load_coverage(coverage_file),
        load_segmentation(segmentation_file),
        load_bedpe(bedpe_file),
    )

# ====== SV cluster helpers ======
def clusters_on_chr(chr_bedpe: pd.DataFrame, chr_name: str) -> dict:
    """
    For a chromosome slice of BEDPE, return:
      { cluster_short: {"pos": [positions on this chr], "n_sv": count of SV rows on this chr} }
      This information is used to draw shaded spans and labels.    
      The for loop goes over every row in chr_bedpe (one structural variant per row) but only looks at the ones that belong to a cluster (because of dropna(subset=["cluster_short"])). 
      For each row:
        It checks which cluster ID (cluster_short) the row belongs to.
        If this is the first time we see that cluster, it creates a new entry for it in the dictionary out.
        It then checks if either side of the variant (chr1 or chr2) is on the chromosome we’re interested in (chr_name). If so, it saves the corresponding start position into that cluster’s "pos" list.
        Finally, it increases the "n_sv" counter by 1, so we know how many variants were assigned to that cluster.
    So, the loop is essentially grouping variants by their cluster ID, keeping track of:
        1) All relevant positions on the given chromosome.
        2) How many variants belong to that cluster in total.
    """
    out = {}
    for _, sv in chr_bedpe.dropna(subset=["cluster_short"]).iterrows():
        cid = sv["cluster_short"]
        if cid not in out:
            out[cid] = {"pos": [], "n_sv": 0}
        # collect positions present on THIS chromosome
        if sv["chr1"] == chr_name:
            out[cid]["pos"].append(int(sv["start1"]))
        if sv["chr2"] == chr_name:
            out[cid]["pos"].append(int(sv["start2"]))
        out[cid]["n_sv"] += 1  # count the row once
    return out

def cluster_colors(cluster_ids) -> dict:
    """
    Stable color per cluster id using tab20 palette.
    It converts the numeric code into an index for the color palette
    It looks up that color from the palette (cmap(...)).
    Stores the result in a dictionary, using the cluster ID (cid) as the key and the color as the value.
    Returns the dictionary mapping each cluster ID to its assigned color.
    """
    cmap = plt.get_cmap("tab20")
    keys = pd.Index(cluster_ids)
    codes = pd.factorize(list(cluster_ids))[0]
    return {cid: cmap(int(code) % 20) for cid, code in zip(keys, codes)}

# ====== plotting helpers ======
def plot_copy_number(ax_cn, chr_cov, chr_seg):
    """
    Plots copy number variation (CNV) data for one chromosome:
        It shows per-bin copy number estimates as a scatter plot.
        It overlays segmentation results (longer regions with uniform copy number) as red lines.
        It adjusts axis labels and grid.
        It returns the maximum x-coordinate (chromosome position) that was drawn, so the caller knows the plotting range.
    """
    # chr_cov contains coverage bins → each row has start, end, and total_cn (the estimated copy number for that bin).
    # If it’s not empty, the code computes the midpoint for each bin
    # This gives an x-coordinate to place each scatter point.

    if not chr_cov.empty:
        chr_cov["mid"] = (chr_cov["start"] + chr_cov["end"]) / 2.0
        ax_cn.scatter(
            chr_cov["mid"], # x-coordinates (midpoints of bins)
            chr_cov["total_cn"], # y-coordinates (copy number estimates)
            s=1, # point size
            alpha=0.6, # transparency of the points
            color="lightblue", 
            label="Coverage (bins)", 
            rasterized=True # helps with large number of points
        )

    # chr_seg contains segments → regions where copy number is constant across multiple bins.
    # For each segment:
    # Draw a horizontal red line at the segment’s copy number value, from its start to its end.
    # Draw two vertical dashed red lines at the start and end positions.
    # These serve as markers showing where the segment boundaries are.

    for _, seg in chr_seg.iterrows(): 
        ax_cn.plot([seg["start"], seg["end"]], # x-coordinates (start to end of segment)
                    [seg["total_cn"], seg["total_cn"]], # y-coordinates (constant copy number)
                   "r-", linewidth=2, alpha=0.8)
        ax_cn.axvline(seg["start"], color="red", linestyle="--", alpha=0.25, linewidth=0.6)
        ax_cn.axvline(seg["end"],   color="red", linestyle="--", alpha=0.25, linewidth=0.6)

    ax_cn.set_ylim(*CN_YLIM) # fixes the y-axis range to a predefined constant CN_YLIM
    ax_cn.set_ylabel("Copy number", labelpad=6)
    ax_cn.grid(True, alpha=0.3)

    # Return the maximum x-coordinate (chromosome position) that was drawn, so the caller knows the plotting range.
    return max(
        chr_cov["end"].max() if not chr_cov.empty else 0,
        chr_seg["end"].max() if not chr_seg.empty else 0
    )

def plot_cluster_spans(ax_sv, ax_cn, clusters_dict, min_svs_per_cluster=2):
    """
    The function plot_cluster_spans highlights structural variant (SV) clusters on two panels (SV panel and copy number panel).
    It:
    - Filters clusters that have enough SVs (at least min_svs_per_cluster).
    - Assigns each cluster a stable color.
    - Shades a rectangular background span (from min position to max position of the cluster) on both panels.
    - Prepares label information (positions + text) so that cluster labels can be added later.
    - Returns both the label records and the color map.
    """
    # Loops through clusters_dict and only keeps clusters where the number of SVs (n_sv) is at least the threshold. This prevents tiny clusters (n = 1) from being shaded or labeled.
    keep = {cid: rec for cid, rec in clusters_dict.items() if rec["n_sv"] >= min_svs_per_cluster } # Looks like this keep = {"A": {"pos": [100, 200, 150], "n_sv": 3}, "B": {"pos": [500, 550], "n_sv": 2}}
    color_map = cluster_colors(list(keep.keys())) # Uses the helper function (cluster_colors) to assign a unique, stable color to each remaining cluster. color_map = { "A": (0.12, 0.47, 0.71, 1.0),  # blue-ish RGBA "B": (1.0, 0.5, 0.05, 1.0) # orange-ish RGBA }

    label_recs = [] # to store label tuples with info for labeling later (xmin, xmax, mid, cluster ID, label text).
    for cid, rec in keep.items(): # Rec is a dict with "pos" (list of positions on this chr) and "n_sv" (count of SVs on this chr).
        xs = rec["pos"] 
        if not xs:  # guard (shouldn't happen)
            continue
        xmin, xmax = min(xs), max(xs) # determine the leftmost and rightmost positions of the cluster on this chromosome.
        # shade on both panels
        c = color_map[cid]
        ax_sv.axvspan(xmin, xmax, color=c, alpha=CLUSTER_ALPHA_SV, zorder=0) # axvspan draws a vertical span (rectangle) between xmin and xmax.
        ax_cn.axvspan(xmin, xmax, color=c, alpha=CLUSTER_ALPHA_CN, zorder=0)
        # prepare label record
        mid = (xmin + xmax) / 2.0
        label_recs.append((xmin, xmax, mid, cid, f"Cluster {cid} (n={rec['n_sv']})"))

    return label_recs, color_map
    # Load bedpe data (structural variants)
    bedpe = (
        pd.read_csv(
            bedpe_file,
            sep="\t",
            header=None,
            dtype=bedpe_dtype
        )
        .rename(columns=str)
        .rename(columns=bedpe_standard_names)
        # Make a short cluster label by taking the part after the last : (so …:10 → "10")
        .assign(
            cluster_short=lambda df: df["cluster_id"]
            .fillna(".")
            .apply(lambda x: x.split(":")[-1] if x != "NA" else x)
        )
        # Sort by chr position
        .sort_values(by=["chr1", "start1", "chr2", "start2"])
        .reset_index(drop=True)
    )
    # 1) cluster_size via value_counts
    counts = bedpe["cluster_short"].value_counts()
    # mapto bedpe
    bedpe["cluster_size"] = bedpe["cluster_short"].map(counts).fillna(0).astype(int)
    # 2) default to "simple"
    bedpe["sv_type"] = "simple"
    # 3) complex if complex_type starts with complex
    bedpe.loc[bedpe["complex_type"].str.startswith("complex", na=False), "sv_type"] = "complex"
    # 4) clustered if not complex and cluster_size >= 2
    bedpe.loc[(bedpe["sv_type"] != "complex") & (bedpe["cluster_size"] >= 2), "sv_type"] = "clustered"

    # Print a quick check: count how many rows have a non-missing cluster_id, and how many have complex_type starting with complex
    print(f"Loaded {bedpe.shape[0]} structural variants, {bedpe['cluster_short'].notna().sum()} with cluster_id, {bedpe['complex_type'].str.startswith('complex').sum()} complex")
    # Print counts per sv_type
    print(bedpe["sv_type"].value_counts())
    # Print counts per complex_type
    print(bedpe["complex_type"].value_counts())
    print(bedpe.head(5))

    return coverage, segmentation, bedpe

    """
    Return DataFrame with columns: chr, pos (TSS), label for the given symbols.
    """
    rel = _RELEASE.get(assembly, 109)
    ens = EnsemblRelease(rel)
    ens.index()  # ensure data is available (downloads/cache on first run)

    rows = []
    for sym in gene_symbols:
        try:
            genes = ens.genes_by_name(sym)
            if not genes:
                continue
            # prefer protein-coding if multiple
            gene = sorted(genes, key=lambda g: g.biotype != "protein_coding")[0]
            chrom = str(gene.contig)
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            # TSS: start on +, end on -
            tss = gene.start if gene.strand == "+" else gene.end
            rows.append({"chr": chrom, "pos": int(tss), "label": sym})
        except Exception as e:
            print(f"[genes] skip {sym}: {e}")
    return pd.DataFrame(rows)

def place_cluster_labels(ax_sv, label_recs):
    """
    Place cluster labels near the top of the SV band using a few 'lanes'
    to avoid overlap.
    Args:
        ax_sv (Axes): Matplotlib Axes for the SV panel
        label_recs (list): List of tuples (xmin, xmax, mid, cid, text) for each cluster to label
    """
    if not label_recs:
        return
    label_recs.sort(key=lambda t: t[0])  # sort clusters by xmin, which is also the genomic position
    lanes_y = [0.95, 0.90, 0.85, 0.80] # relative y positions (0.0–1.0 in axis coordinates) for 4 rows of labels near the top.
    lane_end = [-np.inf] * len(lanes_y) # keep track of the rightmost x position used in each lane to avoid overlap.
    x0, x1 = ax_sv.get_xlim() # get the current x-axis limits of the SV panel.
    pad = (x1 - x0) * 0.02 # 2% of the x-axis width as padding to avoid labels being too close.

    for xmin, xmax, mid, cid, text in label_recs:
        # pick a lane that is free enough; otherwise last lane
        lane = len(lanes_y) - 1
        for i in range(len(lanes_y)):
            if xmin > lane_end[i] + pad: # If the current label’s left edge (xmin) is past that right edge plus a little pad, the lane is “free enough,” so choose it and stop searching.
                lane = i
                break
        y = lanes_y[lane]
        ax_sv.text(
            mid, y, text, ha="center", va="top", fontsize=8, color="black",
            bbox=dict(facecolor="white", alpha=0.65, boxstyle="round,pad=0.2"),
            zorder=3
        )
        lane_end[lane] = max(lane_end[lane], xmax)

    
def plot_svs(ax_sv, ax_cn, chr_bed, chr_name, cluster_color_all):
    """
    Intra-chromosomal SVs (both breakpoints on chr_name): draw a semicircular arc between the two breakpoints and markers at each breakpoint.
    Inter-chromosomal SVs (one breakpoint on chr_name, the partner on some other chromosome): 
        draw a local marker and vertical guide line(s) at the in-chrom breakpoint, and collect partner labels (like chr3:12,345) to place outside the axes, above the panel, stacked to avoid collisions.
    """
    ax_sv.set_ylim(0, 1)
    ax_sv.set_yticks([])
    ax_sv.set_ylabel("SVs", labelpad=10)

    arch_y, arch_h = 0.10, 0.90
    inter_labels = []  # collect (x, "chr:pos") for out-of-axes placement later

    for _, sv in chr_bed.iterrows():
        # color by cluster_short if available; darkblue fallback for singletons
        sv_color = cluster_color_all.get(sv.get("cluster_short"), "darkblue")
        # category and marker from complex_type (col 18)
        # based on a derived SV category from complex_type (e.g., INS/DEL/INV/…); SV_MARKER maps category → matplotlib marker string.
        cat = _sv_category(sv.get("complex_type"))
        marker = SV_MARKER[cat]

        intra = (sv["chr1"] == chr_name) and (sv["chr2"] == chr_name)
        if intra:
            x1, x2 = int(sv["start1"]), int(sv["start2"])
            if x1 == x2:
                continue
            center, width = (x1 + x2) / 2.0, abs(x2 - x1)
            # arc + breakpoint markers
            arc = patches.Arc((center, arch_y), width, arch_h, angle=0, theta1=0, theta2=180,
                              color=sv_color, linewidth=1.5, alpha=0.8, zorder=2)
            ax_sv.add_patch(arc)
            ax_sv.plot(x1, arch_y, marker=marker, markersize=SV_MARKER_SIZE,
                       alpha=0.9, color=sv_color, linestyle="None", zorder=3)
            ax_sv.plot(x2, arch_y, marker=marker, markersize=SV_MARKER_SIZE,
                       alpha=0.9, color=sv_color, linestyle="None", zorder=3)
        else:
            # inter-chrom: draw local stub(s), vertical guides on both panels,
            # and remember partner labels to place outside the axes.
            if sv["chr1"] == chr_name:
                x_local = int(sv["start1"])
                partner = f"{sv['chr2']}:{int(sv['start2']):,}"
                ax_sv.plot(x_local, arch_y, marker=marker, markersize=SV_MARKER_SIZE,
                           alpha=0.9, color=sv_color, linestyle="None", zorder=2)
                ax_sv.axvline(x_local, color=INTER_COLOR, alpha=0.7, linewidth=1.2, zorder=1)
                ax_cn.axvline(x_local, color=INTER_COLOR, alpha=0.4, linewidth=1.0, zorder=0)
                inter_labels.append((x_local, partner))
            if sv["chr2"] == chr_name:
                x_local = int(sv["start2"])
                partner = f"{sv['chr1']}:{int(sv['start1']):,}"
                ax_sv.plot(x_local, arch_y, marker=marker, markersize=SV_MARKER_SIZE,
                           alpha=0.9, color=sv_color, linestyle="None", zorder=2)
                ax_sv.axvline(x_local, color=INTER_COLOR, alpha=0.7, linewidth=1.2, zorder=1)
                ax_cn.axvline(x_local, color=INTER_COLOR, alpha=0.4, linewidth=1.0, zorder=0)
                inter_labels.append((x_local, partner))

    # place partner labels OUTSIDE, *above* the SV axes, stacked per x
    # group all labels by their integer x so labels that share the same genomic x position stack at the same horizontal location.
    if inter_labels:
        grouped = defaultdict(list)
        for x, txt in inter_labels:
            grouped[int(x)].append(txt)
        trans = mtransforms.blended_transform_factory(ax_sv.transData, ax_sv.transAxes)
        y0, dy = 1.02, 0.06  # just above the top spine; stack upward
        for x in sorted(grouped):
            for i, txt in enumerate(grouped[x]):
                ax_sv.text(
                    x, y0 + i*dy, txt, rotation=90, ha="center", va="bottom",
                    color=INTER_COLOR, fontsize=7,
                    bbox=dict(facecolor="white", alpha=0.6, boxstyle="round,pad=0.15"),
                    transform=trans, clip_on=False, zorder=5
                )

def add_legends(ax_cn, ax_sv):
    """CN legend (bins/segments) and SV marker legend (shape → class)."""
    # CN legend
    bin_proxy = Line2D([0], [0], marker="o", linestyle="None", markersize=4,
                       color="lightblue", label="Coverage (bins)")
    seg_proxy = Line2D([0], [0], color="red", lw=2, label="Segments")
    ax_cn.legend(handles=[bin_proxy, seg_proxy], loc="upper right", frameon=False)

    # SV marker legend (neutral black edges so it doesn't imply cluster color)
    shape_handles = [
        Line2D([0],[0], marker=SV_MARKER['td'],    linestyle='None', markersize=SV_MARKER_SIZE,
               markerfacecolor='none', markeredgecolor='black', label='TD/DUP'),
        Line2D([0],[0], marker=SV_MARKER['del'],   linestyle='None', markersize=SV_MARKER_SIZE,
               markerfacecolor='none', markeredgecolor='black', label='DEL'),
        Line2D([0],[0], marker=SV_MARKER['inv'],   linestyle='None', markersize=SV_MARKER_SIZE,
               markerfacecolor='none', markeredgecolor='black', label='INV'),
        Line2D([0],[0], marker=SV_MARKER['inter'], linestyle='None', markersize=SV_MARKER_SIZE,
               markerfacecolor='none', markeredgecolor='black', label='INTER'),
        Line2D([0],[0], marker=SV_MARKER['complex'], linestyle='None', markersize=SV_MARKER_SIZE,
               markerfacecolor='none', markeredgecolor='black', label='Complex'),
    ]
    leg_markers = ax_sv.legend(handles=shape_handles, title="SV type", loc='upper left', frameon=False)
    ax_sv.add_artist(leg_markers)

# ====== per-chromosome orchestrator ======
def plot_chromosome(chr_name, coverage, segmentation, bedpe, ax_sv, ax_cn):
    """
    Inputs:
      - chr_name: e.g. "chr1"
      - coverage/segmentation/bedpe: full DataFrames
      - ax_sv:  matplotlib Axes for top (SVs)
      - ax_cn:  matplotlib Axes for bottom (Copy number)

    Steps:
      1) slice data to this chromosome
      2) plot CN + segments; collect max x
      3) compute shared x-lims including SV coordinates
      4) cluster spans on both panels + labels
      5) SV arcs + markers + inter-chrom guides + partner labels
      6) legends + cosmetics
    """
    chr_cov = coverage[coverage["chr"] == chr_name].copy()
    chr_seg = segmentation[segmentation["chr"] == chr_name].copy()
    chr_bed = bedpe[(bedpe["chr1"] == chr_name) | (bedpe["chr2"] == chr_name)].copy()
    if chr_cov.empty and chr_seg.empty and chr_bed.empty:
        return

    # (1) CN
    cn_max = plot_copy_number(ax_cn, chr_cov, chr_seg)

    # (2) include SV coords in x-limit
    max_pos = cn_max
    if not chr_bed.empty:
        if (chr_bed["chr1"] == chr_name).any():
            max_pos = max(max_pos, int(chr_bed.loc[chr_bed["chr1"] == chr_name, "start1"].max()))
        if (chr_bed["chr2"] == chr_name).any():
            max_pos = max(max_pos, int(chr_bed.loc[chr_bed["chr2"] == chr_name, "start2"].max()))
    ax_cn.set_xlim(0, max_pos)
    ax_sv.set_xlim(0, max_pos)

    # (3) cluster spans + labels
    clust = clusters_on_chr(chr_bed, chr_name)
    label_recs, cluster_color_span = plot_cluster_spans(ax_sv, ax_cn, clust, min_svs_per_cluster=2)
    place_cluster_labels(ax_sv, label_recs)

    # build a color map for ALL clusters on this chr (singletons get darkblue)
    SINGLETON_COLOR = "darkblue"
    cluster_color_all = {cid: cluster_color_span.get(cid, SINGLETON_COLOR) for cid in clust.keys()}

    # (4) SVs
    plot_svs(ax_sv, ax_cn, chr_bed, chr_name, cluster_color_all)

    # (5) cosmetics + legends
    ax_cn.set_xlabel(f"{chr_name} position (bp)")
    ax_cn.ticklabel_format(style="scientific", axis="x", scilimits=(6, 6))
    add_legends(ax_cn, ax_sv)

# ====== main page loop ======
def create_genomic_plot(coverage_file, segmentation_file, bedpe_file,
                        chromosomes=None, output_file=None):
    coverage, segmentation, bedpe = load_data(coverage_file, segmentation_file, bedpe_file)
    if chromosomes is None:
        chromosomes = sort_chromosomes(coverage["chr"].unique())
        default_output = "genomic_plot_All.pdf"
    else:
        # join chromosome names if it’s a list, otherwise just str()
        if isinstance(chromosomes, (list, tuple)):
            chrom_label = "_".join(map(str, chromosomes))
        else:
            chrom_label = str(chromosomes)
        default_output = f"genomic_plot_{chrom_label}.pdf"

    output_file = Path(output_file or default_output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(output_file) as pdf:
        for chr_name in chromosomes:
            fig, (ax_sv, ax_cn) = plt.subplots(
                2, 1, figsize=(15, 8), constrained_layout=False,
                gridspec_kw={"height_ratios": [1.0, 1.0]}
            )
            # leave room above SV for partner labels
            fig.subplots_adjust(top=0.82, hspace=0.35)

            plot_chromosome(chr_name, coverage, segmentation, bedpe, ax_sv, ax_cn)

            fig.suptitle(f"Chromosome {chr_name}", fontsize=16, fontweight="bold")
            pdf.savefig(fig)  # avoid bbox_inches="tight" to not clip outside labels
            plt.close(fig)

    print(f"Plot created! Output saved to: {output_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create aligned copy-number and structural-variant plots."
    )
    parser.add_argument("--coverage", required=True, help="Binned coverage TSV.")
    parser.add_argument("--segmentation", required=True, help="Copy-number segmentation TSV.")
    parser.add_argument("--bedpe", required=True, help="Structural-variant BEDPE file.")
    parser.add_argument(
        "--output",
        default="genomic_plot_All.pdf",
        help="Output PDF path (default: genomic_plot_All.pdf).",
    )
    parser.add_argument(
        "--chromosome",
        dest="chromosomes",
        action="append",
        help="Chromosome to plot; repeat for multiple chromosomes. Defaults to all.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    create_genomic_plot(
        args.coverage,
        args.segmentation,
        args.bedpe,
        chromosomes=args.chromosomes,
        output_file=args.output,
    )


if __name__ == "__main__":
    main()
