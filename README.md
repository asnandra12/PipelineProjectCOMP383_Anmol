# PipelineProjectCOMP383_Anmol

This is the pipeline project for COMP 383. 
It focuses on track 1 and the differential expression of 4 different transcriptomes.
First we obtained the transcriptomes along with some sample data. We then extracted the cds and quantified the TPM in each sample. Finally we blasted to find other virus genes that contain the most most differentially expressed protein from the extracted cds. 
Along with the python file I have also provided some sample data which is a smaller version of the original donor data. 
With the python script you will also need these packages:
  -Biopython
  -Kallisto
  -Sleuth
  -Pandas
When attempting to run the pipeline, please clone the repository. Then please run the python file "PipelineProject.py", it will copy the other files itself and no other action is needed.
