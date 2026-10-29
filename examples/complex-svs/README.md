# SV + CN Visualization Example

This directory contains the notebook and example data for the structural
variant (SV) and copy-number (CN) visualization shown in the repository's
[front-page README](../../README.md).

The reusable version is [`complex_sv_plot.py`](complex_sv_plot.py). The
notebook contains the same plotting code in an interactive format.

This script provides helper functions to visualize **structural variants (SVs)** and **copy number (CN)** information across chromosomes.  
It is designed to produce two aligned panels:  
- **CN panel (ax_cn)**: scatterplot of coverage bins + CNV segments.  
- **SV panel (ax_sv)**: arcs and markers for SV breakpoints, with optional cluster shading and labels.  

---

## Features

### Copy Number (CN) plotting
- **`plot_copy_number(ax_cn, chr_cov, chr_seg)`**  
  - Plots per-bin copy number estimates as a scatter plot.  
  - Adds segmentation tracks as red horizontal lines with vertical boundaries.  
  - Returns the maximum genomic x-position encountered for axis alignment.  

- **`add_legends(ax_cn, ax_sv)`**  
  - Adds two legends:  
    - CN legend (bins vs. segments) on the CN panel.  
    - SV legend (marker shapes for each SV category) on the SV panel.  

---

### Structural Variant (SV) plotting
- **`clusters_on_chr(chr_bedpe, chr_name)`**  
  - Groups SVs by cluster ID for a given chromosome.  
  - Collects positions and SV counts per cluster.  

- **`cluster_colors(cluster_ids)`**  
  - Assigns stable colors per cluster using matplotlib’s *tab20* palette.  

- **`plot_cluster_spans(ax_sv, ax_cn, clusters_dict, min_svs_per_cluster=2)`**  
  - Shades cluster spans on both SV and CN panels.  
  - Prepares label information for later placement.  
  - Returns `(label_recs, color_map)`.  

- **`place_cluster_labels(ax_sv, label_recs)`**  
  - Places cluster labels above the SV panel, using multiple “lanes” to avoid overlap.  

- **`plot_svs(ax_sv, ax_cn, chr_bed, chr_name, cluster_color_all)`**  
  - Draws intra-chromosomal SVs as arcs with markers.  
  - Draws inter-chromosomal SVs as local stubs + vertical guides.  
  - Places partner chromosome labels outside the SV panel, stacked per breakpoint.  

---

## Dependencies

- Python **3.9+**
- Dependencies listed in the repository [`requirements.txt`](../../requirements.txt)
- [pandas](https://pandas.pydata.org/)  
- [matplotlib](https://matplotlib.org/)

## Running the script

From the repository root, install the dependencies and pass the three input
files to the script:

```bash
python -m pip install -r requirements.txt
python examples/complex-svs/complex_sv_plot.py \
  --coverage path/to/binned_coverage.tsv \
  --segmentation path/to/segmentation.tsv \
  --bedpe path/to/structural_variants.bedpe \
  --output genomic_plot.pdf
```

Add `--chromosome chr7` to restrict the output to one chromosome. The input
column layout is documented in the `load_coverage`, `load_segmentation`, and
`load_bedpe` docstrings in the script.

## Example Output

When run with a chromosome slice of coverage (`chr_cov`), segments (`chr_seg`), and structural variants (`chr_bed`), the script produces two aligned panels:

- **CN panel (top):**
  - Light-blue scatter points showing per-bin coverage (copy number estimates).
  - Red horizontal lines showing copy number segments.
  - Vertical dashed red lines marking segment boundaries.

- **SV panel (bottom):**
  - **Intra-chromosomal SVs:** shown as semicircular arcs with markers at both breakpoints.
  - **Inter-chromosomal SVs:** shown as local markers with vertical guide lines extending across both panels.
  - **Cluster spans:** shaded background regions (optional) with labels above the panel.
  - **Partner chromosome labels:** placed outside, above the SV panel, stacked vertically if multiple labels share the same x-position.

Both panels include legends:
- **CN legend** (upper right): “Coverage (bins)” and “Segments”.
- **SV legend** (upper left): marker shapes mapped to SV types (DEL, DUP/TD, INV, INTER, Complex).

The result is a compact visualization where copy number patterns, structural variants, and their clusters can be inspected together.
