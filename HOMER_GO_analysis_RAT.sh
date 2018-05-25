##### Haven't gotten this to work yet, produces empty files ####


################### HOMER (findGO.pl): Find GO terms from genelist input
### Prefers refseq IDs if possible, but can convert SYMBOLs
## http://homer.salk.edu/homer/motif/
## This findGO.pl

#!/bin/bash
# module add perl-scg; module add MEME; module add homer; module add weblogo
# find *.genes | xargs -n1 qsub -V -cwd -l h_vmem=4G -pe shm 1 ./HOMER_GO_analysis.sh #cluster
# find *.genes | xargs -n1 bash HOMER_GO_analysis_RAT.sh #local, analyze all

name=`basename $1 .genes`
findGO.pl $1 rat $name.HOMER_GO -cpu 4
