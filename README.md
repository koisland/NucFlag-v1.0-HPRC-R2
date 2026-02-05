# HPRC Release 2 - NucFlag v1.0

## Getting Started
Clone repo.
```bash
git clone https://github.com/koisland/NucFlag-HPRC-R2-v1.0
cd NucFlag-HPRC-R2-v1.0
```

Install `snakemake`.
```bash
python -m venv venv
source venv/bin/activate
pip install snakemake
```

## Usage
Run the thing. Installs dependencies with `conda`.
```bash
snakemake -np -j 100 --workflow-profile workflow/profiles/lpc
```
