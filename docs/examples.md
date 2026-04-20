# Examples

## Human mtDNA Dataset

```bash
mitocape --input human.csv --outdir results --species human
````

## Mouse Dataset

```bash
mitocape --input mouse.csv --outdir mouse_results --species mouse
```

## Fly Dataset

```bash
mitocape --input fly.csv --outdir fly_results --species fly
```

## Python Workflow

```python
from mitocape import run_analysis

run_analysis(
    input_file="sample.csv",
    outdir="results",
    species="human"
)
```
