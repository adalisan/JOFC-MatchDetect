# 
#$ -cwd 
#$ -j y 
#$ -S /bin/bash 
#$ -N wiki_JOFC_hypTest

echo "running Wiki 2cond "+$SGE_TASK_ID".Rout"


./R/R CMD BATCH  --no-restore  ~/my_documents/DataFusion/src/wiki_hypTest_2cond.R  "run_wiki_JOFC.Rout"

echo "" 
echo "Done at " `date` 

