#!/bin/bash -l

## Variables
FILE=$1
NAME=`basename $1 .genes`

#SBATCH --job-name $NAME.GO
#SBATCH --output=$NAME.GO.out
#SBATCH --mail-user jchap14@stanford.edu
#SBATCH --mail-type=ALL
# Request run time & memory
#SBATCH --time=1:00:00
#SBATCH --mem=1G
#SBATCH --ntasks-per-node=4
#SBATCH --account=mpsnyder
#SBATCH --nodes=1
#SBATCH --export=ALL

################### HOMER (findGO.pl): Find GO terms from genelist input
### Prefers refseq IDs if possible, but can convert SYMBOLs
## http://homer.salk.edu/homer/motif/
## This findGO.pl

# module add perl-scg; module add MEME; module add homer; module add weblogo
# find *.genes | xargs -n1 sbatch -V -cwd -l h_vmem=4G -pe shm 1 ./HOMER_GO_analysis.sh #cluster
# find *.genes | xargs -n1 bash HOMER_GO_analysis_HUMAN.sh #local, analyze all

## commands
cat $FILE | tr -d '"' > $NAME.temp
findGO.pl $NAME.temp human $NAME.HOMER_GO -cpu 4
rm $NAME.temp


