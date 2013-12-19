# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#$ -N wiki_JOFC_hypTest
#$ -pe openmp 4
#$ -v OMP_NUM_THREADS=4

echo "running Wiki 2cond "+$SGE_TASK_ID".Rout"


~/R_2_15_x64 CMD BATCH  --no-restore  ~/projects/JOFC-MatchDetect/src/wiki_hypTest_2cond.R  "run_wiki_JOFC.Rout"

echo "" 
echo "Done at " `date` 

