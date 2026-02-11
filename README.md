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

To generate output manifest:
```bash
snakemake -np -j 40 finalize_output
```

## Outputs
|file|desc|
|-|-|
|`results/final`|Output symlinked directory of files to upload to HPRC S3 bucket.|
|`results/nucflag`|Output directory of files from NucFlag run.|
|`results/nucflag_{dtype}_{out_dtype}.csv*`|Output manifest files per output datatype with format: `sample_id`, `hap`, `assembly_name`, `location`. There are five output datatypes (`bw`, `bed`, `bed_status`, `bed_qv`, and `plots`) across two datatypes (`hifi` and `ont`).`bw` is all bigWig pileup files, `bed` is NucFlag calls, `bed_status` is NucFlag calls summarized by type, `bed_qv` is estimated QV (Phred-like score) by number of calls, and `plots` are the NucFreq plots with self-identity and mapq added as additional tracks.|
