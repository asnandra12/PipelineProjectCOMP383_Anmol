
#load package
library(sleuth)

#read in the table you made describing samples and kallisto ouput
#assign to variable name stab
stab = read.table("sleuth_table.txt", header=TRUE, stringsAsFactors = FALSE)


#initialize slueth object using sleuth_prep function from sleuth library
so = sleuth_prep(stab)

#fit a model comparing the two conditions
so = sleuth_fit(so,~condition,'full')

#fit the reduced model to compare in the likelihood ratio test
so = sleuth_fit(so,~1,'reduced')

#perform the likelihood ratio test for differential expressions between conditions
so = sleuth_lrt(so,'reduced','full')

#load the dplyr package for data.frame filtering
library(dplyr)

#extract the test results from the sleuth object
sleuth_table = sleuth_results(so,'reduced:full','lrt',show_all = FALSE)

#filter the most significant results(FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table,qval<=0.05) |> dplyr::arrange(pval)

write.table(sleuth_significant,file="fdr_results.txt",quote=FALSE,row.names=FALSE)

