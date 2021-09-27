#!/bin/bash
# DESCRIPTION
#    covid19_wgs_consensus: COVID-19 WGS consensus pipeline.
#    
# IMPLEMENTATION
#    author   S?ren Karst (sorenkarst@gmail.com)
#    license  GNU General Public License

### Description ----------------------------------------------------------------

USAGE="
-- covid19_wgs: COVID-19 WGS consensus and QC pipeline.

usage: covid19_wgs [-h ] ( -i dir -o dir -b string -t integer)
where:
    -h   Show this help text.
    -i   Parent folder containing barcode subfolder with fastq file(s).
    -o   Output folder and working directory.
    -b   Space delimited string og used barcode IDs. IDs should match names of
         barcode subfolders. Bash expansion possible (eg. barcode{25..32}
         barcode{34..42}). If unused barcodes are registered, they will be
         partially processed and flagged in the QC.
    -t   Number of CPU threads.
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzdi:o:b:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    i) DATA_DIR=$OPTARG;;
    o) ANALYSIS_DIR=$OPTARG;;
    b) BARCODE_SUBSET=$(eval echo "$OPTARG" | sed 's/ /|/g');;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${DATA_DIR+x} ]; then echo "-i $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${ANALYSIS_DIR+x} ]; then echo "-o $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${BARCODE_SUBSET+x} ]; then echo "-b $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${THREADS+x} ]; then echo "-t $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${DATA_DIR+x} ]; then echo "-i $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;

# Conda environments dir
export $CONDA_DIR="$(dirname "$(readlink -f "$0")" | sed 's|envs/covid19_wgs|envs/g')"

# Activate conda environment
ARTIC_VERSION=$(artic -v)

# Fastq
mkdir $ANALYSIS_DIR/fastq

find \
  $DATA_DIR/ \
  -maxdepth 1 \
  -mindepth 1 \
  -type d |\
parallel \
  -j $THREADS \
  --rpl '{bid} s:.*/::' \
  "cp -r {} $ANALYSIS_DIR/fastq/"

# Filter
mkdir $ANALYSIS_DIR/filter

find \
  $ANALYSIS_DIR/fastq \
  -type d |\
grep -E "$BARCODE_SUBSET" |\
parallel \
  -j $THREADS \
  --rpl '{bid} s:.*/::' \
  "artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory {} --output $ANALYSIS_DIR/filter/{bid}_filt.fastq"
  
# Consensus
mkdir $ANALYSIS_DIR/consensus

find \
  $ANALYSIS_DIR/filter \
  -type f \
  -name "barcode*filt.fastq" |\
parallel \
  -j $THREADS \
  --rpl '{bid} s:.*/::;s:_filt.fastq::' \
  "mkdir $ANALYSIS_DIR/consensus/{bid};
  cd $ANALYSIS_DIR/consensus/{bid};
  artic minion \\
  --medaka \\
  --medaka-model r941_min_high_g360 \\
  --normalise 100 \\
  --threads 1 \\
  --scheme-directory $CONDA_DIR/artic/artic-ncov2019/primer_schemes \\
  --read-file {} \\
  nCoV-2019/V3 {bid}"

find \
  $ANALYSIS_DIR/consensus \
  -name "*.consensus.fasta" \
  -exec cat {} \; \
  > $ANALYSIS_DIR/consensus.fa

# Coverage

# Repeat align trim without normalization
find \
  $ANALYSIS_DIR/consensus/*/ \
  -maxdepth 1 \
  -type f \
  -regextype posix-extended \
  -regex '.*barcode[0-9][0-9].sorted.bam' | \
parallel \
  -j $THREADS \
  --rpl '{sid} s:.*/::;s:.sorted.*::' \
  "
  align_trim \\
    --normalise 1000000 \\
    $CONDA_DIR/artic/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed \\
    --remove-incorrect-pairs \\
    --report $ANALYSIS_DIR/consensus/{sid}/{sid}.alignreport_nonormalize.txt \\
    < {} 2> \\
    $ANALYSIS_DIR/consensus/{sid}/{sid}.alignreport_nonormalize.er |\\
    samtools sort -T {sid} - -o $ANALYSIS_DIR/consensus/{sid}/{sid}.primertrimmed_nonormalize.rg.sorted.bam
  samtools index $ANALYSIS_DIR/consensus/{sid}/{sid}.primertrimmed_nonormalize.rg.sorted.bam
  "

cat \
  $CONDA_DIR/artic/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.insert.bed \
  > $ANALYSIS_DIR/coverage.bed
  
  
find \
  $ANALYSIS_DIR/consensus/*/ \
  -maxdepth 1 \
  -type f \
  -regextype posix-extended \
  -regex '.*barcode[0-9][0-9].primertrimmed_nonormalize.rg.sorted.bam' | \
parallel \
  -j $THREADS \
  --rpl '{sid} s:.*/::;s:.primertrimmed_nonormalize.rg.sorted.bam::' \
  "bedtools coverage -d -a $ANALYSIS_DIR/coverage.bed -b {} | sed 's/^/{sid}\t/g'" > $ANALYSIS_DIR/coverage.tsv


# Yield
gz_check (){
  if [ "${1##*.}" = "gz" ]; then
    gzip -cd "$1" 
  else
    cat "$1"
  fi
}

export -f gz_check

find \
  $ANALYSIS_DIR/fastq \
  -type f \
  -name "*fastq*" \
  -exec bash -c 'gz_check "$@"' bash {} \; |\
gawk \
  -v OUT="$ANALYSIS_DIR" \
  '
  NR%4==1 {
    DATE=substr($5, 12, 10)
    TIME=substr($5, 23, 5)
    BARCODE=substr($9, 9)
    getline
    LEN=length($0)
    # Create array
    interval[BARCODE"\t"DATE" "TIME]["reads"]++
    interval[BARCODE"\t"DATE" "TIME]["bases"]+=LEN
    len[BARCODE"\t"LEN]++
  }
  END{
    print "barcode_id\tdate_time\treads_n\tyield_bp" > OUT"/yield.tsv"
    for (i in interval){
      print i "\t" interval[i]["reads"] "\t" interval[i]["bases"] > OUT"/yield.tsv"
    }
    print "barcode_id\tlength_bp\tcount" > OUT"/length.tsv"
    for (i in len){
      print i "\t" len[i] > OUT"/length.tsv"
    }
  }
  '
  
# Pangolin classification
conda activate $CONDA_DIR/covid19_wgs_pangolin
PANGOLIN_VERSION=$(pangolin -v)
PANGOLEARN_VERSION=$(pangolin -pv)

pangolin \
  $ANALYSIS_DIR/consensus.fa \
  --outdir $ANALYSIS_DIR \
  --outfile pangolin.csv

# SNP calling and visualization
conda activate $CONDA_DIR/covid19_wgs_nextclade
NXCL_VERSION=$(nextclade --version)
NXCL_REF='MN908947'
NXCL_TAG='2021-06-25T00:00:00Z'
NXCL_DATA_VERSION=$()

mkdir $ANALYSIS_DIR/nextclade

nextclade dataset get \
  --name='sars-cov-2' \
  --reference=$NXCL_REF \
  --tag= $NXCL_TAG\
  --output-dir=$ANALYSIS_DIR/nextclade/sars-cov-2_nxcl_db

nextclade run\
  --input-fasta $ANALYSIS_DIR/consensus.fa \
  --input-dataset $ANALYSIS_DIR/nextclade/sars-cov-2_nxcl_db \
  --output-dir $ANALYSIS_DIR/nextclade \
  --output-tsv $ANALYSIS_DIR/nextclade.tsv \
  -j $THREADS
  
# SNP visualization
conda activate $CONDA_DIR/covid19_wgs_nextstrain
NXST_VERSION=$(nextstrain --version)

# Remove poor sequences
gawk \
  '
    NR%2==1{
      HEAD=$0
      getline
      SEQN=$0
      gsub("[ATCG-]", "", SEQN)
      if(length(SEQN) < 3000){
        print HEAD "\n" $0
      }
    }
  ' \
  $ANALYSIS_DIR/nextclade/consensus.aligned.fasta \
  > $ANALYSIS_DIR/nextclade/consensus.aligned.hgmq.fasta

# Align sequences
augur tree \
  --alignment $ANALYSIS_DIR/nextclade/consensus.aligned.hgmq.fasta \
  --output $ANALYSIS_DIR/consensus_tree.nwk \
  --nthreads $THREADS &> $ANALYSIS_DIR/log.out
  
# Protocol specification file
get_latest_release() {
  curl --silent "https://api.github.com/repos/$1/releases/latest" | # Get latest release from GitHub api
    grep '"tag_name":' |                                            # Get tag line
    sed -E 's/.*"([^"]+)".*/\1/'                                    # Pluck JSON value
}

printf "process\tversion_used\tversion_available
lab_protocol\tNA\tNA
primer_scheme\tArtic V3\tNA
consensus_pipeline\t$ARTIC_VERSION\tartic $(get_latest_release "artic-network/fieldbioinformatics")
pangolin\t$PANGOLIN_VERSION\tpangolin $(get_latest_release "cov-lineages/pangolin")
pangolearn\t$PANGOLEARN_VERSION\tpangoLEARN $(get_latest_release "cov-lineages/pangoLEARN")
nextstrain-cli\t$NXST_VERSION\tNA
nextclade\tNextclade $NXCL_VERSION\tNextclade $(get_latest_release "nextstrain/nextclade")
nextclade_data\tRef:$NXCL_REF;Tag:$NXCL_TAG\tNA
" > $ANALYSIS_DIR/specifications.tsv

