# Plot Types

MitoCAPE can generate multiple visual outputs.

## Available Plot Types

| Plot Name | Description |
|----------|-------------|
| standard | Mutation category barplot |
| top_positions | Most edited positions |
| bubble | Position vs percentage bubble plot |
| pie | Silent vs amino-acid change pie chart |
| donut | Donut version of pie chart |
| histogram | Per-sample distribution |
| line | Per-sample trend plot |
| sequence_logo | Sequence preference logo |
| gene_mut | Gene-level mutation plot |
| gene_codon | Codon-level impact plot |

# How to Generate Plots

Plots are created automatically when running MitoCAPE unless disabled.

## Example 1: Generate Default Plots

```bash
mitocape --input sample.csv --outdir results --species human
```
This typically generates:

standard plot
top edited positions plot
bubble plot

Saved in:

results/
Example 2: Generate Only Bubble Plot
```bash
mitocape --input sample.csv --outdir results --species human --plot-types bubble
```
Output:

results/sample_amino_acid_changes_bubble_plot.pdf
Example 3: Generate Standard + Pie + Donut
```bash
mitocape --input sample.csv --outdir results --species human --plot-types standard pie donut
```
Useful for summaries in presentations or publications.

Python API Example

Generate a bubble plot directly from Python:

from mitocape import run_analysis

```python
run_analysis(
    input_file="sample.csv",
    outdir="results",
    species="human",
    plot_types=["bubble"]
)
```
Sequence Logo Example

To generate sequence context preferences:
```bash
mitocape --input sample.csv --outdir results --species human --plot-types sequence_logo
```
Or in Python:
```python
run_analysis(
    input_file="sample.csv",
    outdir="results",
    species="human",
    plot_types=["sequence_logo"],
    sequence_logo_threshold=10
)
```
This creates a motif/logo using edits above the threshold.

Gene-Level Plot Example
```bash
mitocape --input sample.csv --outdir results --species human --plot-types gene_mut
```

Used to visualise mutation burden across genes.

---

## Output Format

Plots are saved as PDF files by default.

---
Use Multiple Plots Together
```bash
--plot-types standard bubble pie
```
Faster Analysis Without Plots
```bash
--no-plots
```
Focus on Stronger Editing Events
```bash
--value-threshold 5
```
Only edits ≥5% are used in selected plots.
__
Recommended Starter Command
```bash
mitocape --input sample.csv --outdir results --species human --plot-types standard bubble top_positions
```
This gives the most informative first overview.
