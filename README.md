# covid19_wgs
Pipeline for COVID-19 WGS consensus generation and QC. Minimalistic wrapper around ARTIC pipeline.

## Installation

* Install [miniconda3](https://docs.conda.io/en/latest/miniconda.html).
* Install covid19_wgs conda environment.
```bash
# Install env
conda create -n covid19_wgs

# Clone covid19_wgs
conda activate covid19_wgs
git clone https://github.com/SorenKarst/covid19_wgs.git $CONDA_PREFIX/covid19_wgs

# Create link
mkdir -p $CONDA_PREFIX/bin
find \
  $CONDA_PREFIX/covid19_wgs/ \
  -name "*.sh" \
  -exec chmod +x {} \;  
ln -s \
  $CONDA_PREFIX/covid19_wgs/covid19_wgs.sh \
  $CONDA_PREFIX/bin/covid19_wgs
```
