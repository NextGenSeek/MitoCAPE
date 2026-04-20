# Python API

## Import

```python
from mitocape import run_analysis
````

## Basic Example

```python
results = run_analysis(
    input_file="sample.csv",
    outdir="results",
    species="human"
)
```

## Return Object

```python
{
    "processed_file": "...",
    "plots": [...],
    "summary": ...
}
```

## Main Parameters

| Parameter       | Description                  |
| --------------- | ---------------------------- |
| input_file      | Input CSV                    |
| outdir          | Output folder                |
| species         | human / mouse / fly          |
| value_threshold | Minimum percentage threshold |
| make_plots      | True / False                 |
| plot_types      | List of plot names           |

## Example with Custom Plots

```python
run_analysis(
    input_file="sample.csv",
    outdir="results",
    species="mouse",
    plot_types=["standard", "bubble"]
)
```
