# covid19_wgs
Pipeline for COVID-19 WGS consensus generation and QC. Minimalistic wrapper around ARTIC pipeline.

## Installation

* Install [miniconda3](https://docs.conda.io/en/latest/miniconda.html).
* Install [mamba package manager](https://github.com/mamba-org/mamba).
* Create conda environments.
```bash
# NOTES:
# Environments are split to ensure easy update.
# Artic should be updated when major changes in ONT chemistry or basecalling occures.
# Pangolin should be regularly updated for correct and current classification.
# Nextclade should be updated when major bug fixes og changes are implemented.
# Nextstrain probably won't need updating. It is only used for tree visualization.

# Install base environment with artic
mamba create \
  -n covid19_wgs \
  artic=1.2.1 \
  bedtools=2.30.0 \
  parallel=20170422 \
  -c bioconda \
  -c conda-forge
  
# Install nextstrain environment
mamba create \
 -n covid19_wgs_nextstrain \
 nextstrain=20200304 \
 -c bioconda 
  
# Install nextclade environment 
mamba create \
 -n covid19_wgs_nextclade \
 nextclade=1.3.0 \
 -c bioconda
  
# Install pangolin environment
mamba create \
 -n covid19_wgs_pangolin \
 pangolin=3.1.11 \
 -c bioconda
 
# Clone covid19_wgs and artic
conda activate covid19_wgs
git clone \
  https://github.com/SorenKarst/covid19_wgs.git \
  $CONDA_PREFIX/covid19_wgs

git clone \
  https://github.com/artic-network/artic-ncov2019.git \
  $CONDA_PREFIX/artic-ncov2019
  
# Copy custom schemes into artic-ncov2019
cp \
  -r \
  $CONDA_PREFIX/covid19_wgs/schemes/V1M \
  $CONDA_PREFIX/artic-ncov2019/primer_schemes/nCoV-2019/  

# Make scripts excutable
mkdir -p $CONDA_PREFIX/bin
find \
  $CONDA_PREFIX/covid19_wgs/ \
  -name "*.sh" \
  -exec chmod +x {} \;  
ln -s \
  $CONDA_PREFIX/covid19_wgs/covid19_wgs.sh \
  $CONDA_PREFIX/bin/covid19_wgs
```
