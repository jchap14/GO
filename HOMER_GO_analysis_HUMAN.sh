#!/bin/bash
##### Trim FASTQs for quality & adaptors. FASTQC output generated also

##### for x in `/bin/ls *.untrimmed.R1.fq.gz` ; do bash bash HOMER_GO_analysis_HUMAN.sh $x; done

## add modules (in Conda env, so modules not required)
# module add perl-scg; module add MEME; module add homer; module add weblogo

## define variables
FILE=$1
NAME=`basename $1 .genes`

## write a tempscript to be looped over
cat > $NAME.tempscript.sh << EOF
#!/bin/bash -l
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

## commands
cat $FILE | tr -d '"' > $NAME.temp
findGO.pl $NAME.temp human $NAME.HOMER_GO -cpu 4
rm $NAME.temp
##
EOF

## qsub then remove the tempscript
sbatch $NAME.tempscript.sh 
sleep 1
rm $NAME.tempscript.sh

