#!/bin/sh

set -e

log() {
        tput bold
        tput setaf 2
        printf "%12s " $1
        tput sgr0
        shift 1
        echo $@
}

error() {
        tput bold
        tput setaf 1
        printf "%12s " $1
        tput sgr0
        shift 1
        echo $@
}


if [ ! -e "$1" ]; then
  error Failed to locate input index file
  exit 1
elif [ ! -e "$2" ]; then
  error Failed to locate input HMM file
  exit 1
fi


JOB_ID=$(uuidgen)
LOCALDIR=$(dirname $(realpath $0))
INPUTDIR=$(dirname $(dirname $(dirname $(realpath $0))))
HMM_PATH="$2"
HMM_NAME=$(basename $HMM_PATH .hmm.gz)
FAIDX="$1"
FA=$(echo $FAIDX | sed 's/.fai//g')

OUT_FILE=$(echo $1 | sed 's/.fna.fai//g').$HMM_NAME.features.tsv
log Result file will be created at $OUT_FILE

# create temporary folder for list files
mkdir -p $LOCALDIR/tmp
mkdir -p $LOCALDIR/tmp/wheels

# Load modules
log Loading required modules
module load Python/3.8.6-GCCcore-10.2.0

# Find every FASTA file in the data directory, and cache that for analysis
log Collecting all sequences in input directory
cut -f1 -d $'\t' "$1" > $LOCALDIR/tmp/${JOB_ID}.files.txt
FA_COUNT=$( cat $LOCALDIR/tmp/${JOB_ID}.files.txt | wc -l )
log Found $FA_COUNT contigs to process

# Deduce the job count and bin size
#JOB_COUNT=$( python -c "print(min(2000, $FA_COUNT))" )
#BIN_SIZE=$( python -c "import math; print(math.ceil(max($FA_COUNT/$JOB_COUNT, 1)))" )
#JOB_COUNT=$( python -c "import math; print(math.ceil($FA_COUNT / $BIN_SIZE))" )
BIN_SIZE=1
JOB_COUNT=$FA_COUNT
log Using $JOB_COUNT jobs

# Download wheels
pip wheel https://github.com/zellerlab/GECCO/archive/dev.zip -w $LOCALDIR/tmp/wheels

# ---

mkdir -p $LOCALDIR/out/${JOB_ID}
mkdir -p $LOCALDIR/logs/out/$JOB_ID
mkdir -p $LOCALDIR/logs/err/$JOB_ID

# ---
JOB_SCRIPT="$LOCALDIR/tmp/${JOB_ID}.job.sh"
log Writig job script to \`$JOB_SCRIPT\`
echo '#!/bin/sh'                                                              >  "$JOB_SCRIPT"
echo '#SBATCH -A zeller'                                                      >> "$JOB_SCRIPT"
echo '#SBATCH -t 02:00:00'                                                    >> "$JOB_SCRIPT"
echo '#SBATCH -n 1'                                                           >> "$JOB_SCRIPT"
echo '#SBATCH --mail-type ALL'                                                >> "$JOB_SCRIPT"
echo '#SBATCH --mem 2G'                                                       >> "$JOB_SCRIPT"
echo '#SBATCH --oversubscribe'                                                >> "$JOB_SCRIPT"
echo "#SBATCH -o $LOCALDIR/logs/out/$JOB_ID/%a.out"                           >> "$JOB_SCRIPT"
echo "#SBATCH -e $LOCALDIR/logs/err/$JOB_ID/%a.err"                           >> "$JOB_SCRIPT"
echo "#SBATCH --gres=tmp:512M"                                                >> "$JOB_SCRIPT"
echo "#SBATCH --tmp=512M"                                                     >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo 'module load Python/3.8.6-GCCcore-10.2.0 SAMtools'                       >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo 'set -e'                                                                 >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo "FA=$FA"                                                                 >> "$JOB_SCRIPT"
echo "HMM_NAME=$HMM_NAME"                                                     >> "$JOB_SCRIPT"
echo "FILELIST=$LOCALDIR/tmp/$JOB_ID.files.txt"                               >> "$JOB_SCRIPT"
echo "OUTDIR=$LOCALDIR/out/$JOB_ID"                                           >> "$JOB_SCRIPT"
echo "INPUTDIR=$INPUTDIR"                                                     >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo "BIN_SIZE=$BIN_SIZE"                                                     >> "$JOB_SCRIPT"
echo 'OFFSET=$(( $BIN_SIZE * ($SLURM_ARRAY_TASK_ID - 1) ))'                   >> "$JOB_SCRIPT"
echo 'BIN_SEQS=$(tail $FILELIST -n +$OFFSET | head -n $BIN_SIZE)'             >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo 'virtualenv $TMPDIR/venv'                                                >> "$JOB_SCRIPT"
echo 'source $TMPDIR/venv/bin/activate'                                       >> "$JOB_SCRIPT"
echo "pip install --no-cache-dir $LOCALDIR/tmp/wheels/*"                      >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo "cp -v $HMM_PATH \$TMPDIR/${HMM_NAME}.hmm.gz"                            >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo 'for seq in $BIN_SEQS; do'                                               >> "$JOB_SCRIPT"
echo '  base=$(basename $seq)'                                                >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo '	echo Copying $base sequence to $TMPDIR'                               >> "$JOB_SCRIPT"
echo '	samtools faidx $FA $seq -o $TMPDIR/$base'                             >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo '  echo Running GECCO on $base'                                          >> "$JOB_SCRIPT"
echo '  mkdir -p "$TMPDIR/out/"'                                              >> "$JOB_SCRIPT"
echo '  gecco -vv annotate -j1 --e-filter 1000 --genome $TMPDIR/$base --hmm $TMPDIR/${HMM_NAME}.hmm.gz -o $TMPDIR/out/$base.$HMM_NAME.features.tsv' >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo '	echo Copying $base results to $OUTDIR'                                >> "$JOB_SCRIPT"
echo '  mkdir -p "$OUTDIR/$SLURM_ARRAY_TASK_ID"'                              >> "$JOB_SCRIPT"
echo '  cp $TMPDIR/out/$base.$HMM_NAME.features.tsv -t "$OUTDIR/$SLURM_ARRAY_TASK_ID/"'        >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"
echo 'done'                                                                   >> "$JOB_SCRIPT"
echo ''                                                                       >> "$JOB_SCRIPT"


# ---

log Submitting jobs to the cluster
sbatch -W --array=1-$JOB_COUNT $JOB_SCRIPT
#sbatch -W --qos=high --array=1 $JOB_SCRIPT

head -n1 $(ls -1 $LOCALDIR/out/$JOB_ID/*/*.$HMM_NAME.features.tsv | head -n1) > $OUT_FILE
for i in $(seq 1 $JOB_COUNT); do
  file=$(ls $LOCALDIR/out/$JOB_ID/$i/*.$HMM_NAME.features.tsv)
  log Collecting results from $(realpath $file --relative-to $LOCALDIR/out/$JOB_ID)
  tail -n+2 $file >> $OUT_FILE
done

log Finished creating table $OUT_FILE
log Removing temporary files
rm -rd $LOCALDIR/out/$JOB_ID
