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
git clone \
  https://github.com/cov-lineages/pangolin.git\
  --branch v3.1.14
cd pangolin
mamba env create -f environment.yml -n covid19_wgs_pangolin
conda activate covid19_wgs_pangolin
pip install .
 
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
  $CONDA_PREFIX/covid19_wgs/custom_schemes/. \
  $CONDA_PREFIX/artic-ncov2019/primer_schemes/nCoV-2019/ 
  
# Fix naming error in artic V4 scheme
for FILE in $CONDA_PREFIX/artic-ncov2019/primer_schemes/nCoV-2019/V4/SARS-CoV-2*
do
  FILE_FORMAT=$(echo "$FILE" | sed 's|SARS-CoV-2|nCoV-2019|')
  mv "$FILE" "$FILE_FORMAT"
done


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
