#!/bin/bash
# DESCRIPTION
#    covid19_wgs: COVID-19 WGS consensus pipeline.
#    
# IMPLEMENTATION
#    author   Soren Karst (sorenkarst@gmail.com)
#    license  GNU General Public License

# Pre-run setup

### Conda
CONDA_DIR=$(dirname $(readlink -f "$0") | sed 's|envs/covid19_wgs/covid19_wgs|envs|g')
CONDA=$(echo $CONDA_DIR | sed 's|envs$|etc/profile.d/conda.sh|g')
source $CONDA

### Primer schemes
SCHEME_LIST=$(
  find \
    $CONDA_DIR/covid19_wgs/artic-ncov2019/primer_schemes/nCoV-2019/ \
    -maxdepth 1 \
    -mindepth 1 \
    -type d |\
  sed 's|.*/||'
)

### Description ----------------------------------------------------------------

USAGE="
-- covid19_wgs: COVID-19 WGS consensus and QC pipeline for ONT data.

usage: covid19_wgs [-h ] ( -i dir -o dir -b string - p string -s string -t integer)
where:
    -h   Show this help text.
    -i   Parent folder containing barcode subfolder with fastq file(s).
    -o   Output folder and working directory.
    -b   Space delimited string of used barcode IDs. IDs should match names of
         barcode subfolders. Bash expansion possible (eg. \"barcode{25..32} barcode{34..42}\").
         If data from unused barcodes are registered, they will be partially
         processed and flagged in the QC.
    -p   Library protocol {artic_V3, midnight_V1M, ...}.
         artic_V3: intact amplicons, expect raw reads to be 400-700 bp.
         midnight_V1M: fragmented amplicons, expect raw reads to be 250-1500 bp.
         ...: Custom library protocol in the format \"primer_scheme;length_min;length_max\"
         eg. \"V4;250;400\"
    -a   Optional arguments to consensus pipeline.
         ONT default: \"--medaka --medaka-model r941_min_high_g360\"
         ILM compatible: \"--bwa --skip-nanopolish --no-longshot\" (Only use for testing)
    -t   Number of CPU threads.
    
Available primer schemes: $(echo $SCHEME_LIST)
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzdi:o:b:p:at:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    i) DATA_DIR=$OPTARG;;
    o) ANALYSIS_DIR=$OPTARG;;
    b) BARCODE_SUBSET=$(eval echo "$OPTARG" | sed 's/ /|/g');;
    p) PROTOCOL=$OPTARG;;
    a) ARTIC_ARGS=$OPTARG;;
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
if [ -z ${PROTOCOL+x} ]; then echo "-p $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${ARTIC_ARGS+x} ]; then ARTIC_ARGS="--medaka --medaka-model r941_min_high_g360"; fi;
if [ -z ${THREADS+x} ]; then echo "-t $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;


# Setup protocol settings
case $PROTOCOL in
  artic_V3)
    LEN_MIN=400
    LEN_MAX=700
    PRIMER_SCHEME=V3
    ;;
  midnight_V1M)
    LEN_MIN=250
    LEN_MAX=1500
    PRIMER_SCHEME=V1M
    ;;
  *\;*\;*)
    LEN_MIN=$(echo "$PROTOCOL" | cut -d";" -f2)
    LEN_MAX=$(echo "$PROTOCOL" | cut -d";" -f3)
    PRIMER_SCHEME=$(echo "$PROTOCOL" | cut -d";" -f1)
  ;;
  *) printf "Invalid format for -p: $PROTOCOL\n" >&2
  echo ""
  echo "$USAGE"
  exit 1
  ;;
esac

# Setup analysis folder
PWD=$(pwd)
mkdir $ANALYSIS_DIR

# Artic version
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
  "artic guppyplex \\
  --skip-quality-check \\
  --min-length $LEN_MIN \\
  --max-length $LEN_MAX \\
  --directory {} \\
  --output $ANALYSIS_DIR/filter/{bid}_filt.fastq"
  
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
  --medaka --medaka-model r941_min_high_g360 \\
  --normalise 100 \\
  --threads 1 \\
  --scheme-directory $CONDA_DIR/covid19_wgs/artic-ncov2019/primer_schemes \\
  --read-file $PWD/{} \\
  nCoV-2019/$PRIMER_SCHEME {bid}"

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
    --normalise 100000 \\
    $CONDA_DIR/covid19_wgs/artic-ncov2019/primer_schemes/nCoV-2019/$PRIMER_SCHEME/nCoV-2019.scheme.bed \\
    --remove-incorrect-pairs \\
    --report $ANALYSIS_DIR/consensus/{sid}/{sid}.alignreport_nonormalize.txt \\
    < {} 2> \\
    $ANALYSIS_DIR/consensus/{sid}/{sid}.alignreport_nonormalize.er |\\
    samtools sort -T {sid} - -o $ANALYSIS_DIR/consensus/{sid}/{sid}.primertrimmed_nonormalize.rg.sorted.bam
  samtools index $ANALYSIS_DIR/consensus/{sid}/{sid}.primertrimmed_nonormalize.rg.sorted.bam
  "

cat \
  $CONDA_DIR/covid19_wgs/artic-ncov2019/primer_schemes/nCoV-2019/$PRIMER_SCHEME/nCoV-2019.insert.bed \
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
    FLOWCELL=substr($6, 14)
    getline
    LEN=length($0)
    # Create array
    interval[BARCODE"\t"DATE" "TIME"\t"FLOWCELL]["reads"]++
    interval[BARCODE"\t"DATE" "TIME"\t"FLOWCELL]["bases"]+=LEN
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
      if(length(SEQN) < 7000){
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

# Copy QC report template
cp \
  $CONDA_DIR/covid19_wgs/covid19_wgs/covid19_wgs_qc.rmd \
  $ANALYSIS_DIR/

# Protocol specification file
get_latest_release() {
  curl --silent "https://api.github.com/repos/$1/releases/latest" | # Get latest release from GitHub api
    grep '"tag_name":' |                                            # Get tag line
    sed -E 's/.*"([^"]+)".*/\1/'                                    # Pluck JSON value
}

# Print specifications

if [ "$PROTOCOL" = "$PRIMER_SCHEME;$LEN_MIN;$LEN_MAX" ]; then
  PRINT_PROTOCOL="$PROTOCOL"
else
  PRINT_PROTOCOL="$PROTOCOL ($PRIMER_SCHEME;$LEN_MIN;$LEN_MAX)"
fi
printf "process\tversion_used\tversion_available
sequencing_run_id\t\tNA
bioinformatics_protocol\t$PRINT_PROTOCOL\tNA
consensus_pipeline\t$ARTIC_VERSION\tartic $(get_latest_release "artic-network/fieldbioinformatics")
consensus_arguments\t$ARTIC_ARGS\tNA
pangolin\t$PANGOLIN_VERSION\tpangolin $(get_latest_release "cov-lineages/pangolin")
pangolearn\t$PANGOLEARN_VERSION\tpangoLEARN $(get_latest_release "cov-lineages/pangoLEARN")
nextstrain-cli\t$NXST_VERSION\tNA
nextclade\tNextclade $NXCL_VERSION\tNextclade $(get_latest_release "nextstrain/nextclade")
nextclade_data\tRef:$NXCL_REF;Tag:$NXCL_TAG\tNA
" > $ANALYSIS_DIR/specifications.tsv

