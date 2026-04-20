# Command Line Usage

## Basic Run

```bash
mitocape --input sample.csv --outdir results --species human
````

## Required Arguments

| Argument | Description    |
| -------- | -------------- |
| --input  | Input CSV file |
| --outdir | Output folder  |

## Optional Arguments

| Argument          | Description                  |
| ----------------- | ---------------------------- |
| --species         | human / mouse / fly          |
| --plot-types      | Select plots                 |
| --value-threshold | Minimum percentage threshold |
| --fasta-path      | Custom FASTA                 |
| --gene-map-path   | Custom gene annotation CSV   |

## Example: Fly Analysis

```bash
mitocape --input fly.csv --outdir fly_results --species fly
```

## Example: Custom Threshold

```bash
mitocape --input sample.csv --outdir results --species human --value-threshold 5
```

## Example: Specific Plots

```bash
mitocape --input sample.csv --outdir results --plot-types standard bubble top_positions
```
