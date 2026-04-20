# Installation

## Install from GitHub

```bash
pip install git+https://github.com/NextGenSeek/MitoCAPE.git
````

## Recommended: Virtual Environment

```bash
python -m venv .venv
source .venv/bin/activate      # Mac/Linux
.venv\Scripts\activate         # Windows
pip install git+https://github.com/NextGenSeek/MitoCAPE.git
```

## Test Installation

```bash
mitocape --help
```

or

```python
from mitocape import run_analysis
```

## Dependencies

Installed automatically:

* pandas
* numpy
* matplotlib
* seaborn
* biopython
* logomaker

## Troubleshooting

### Command not found

Use:

```bash
python -m mitocape.cli --help
```

### Upgrade package

```bash
pip install --upgrade git+https://github.com/NextGenSeek/MitoCAPE.git
```
