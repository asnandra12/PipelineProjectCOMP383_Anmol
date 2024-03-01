# PipelineProjectCOMP383_Anmol

This is the pipeline project for Dr. Wheeler's COMP 383. 
It focuses on track 1 and the differential expression of 4 different transcriptomes.
First we obtained the transcriptomes along with some sample data. We then extracted the cds and quantified the TPM in each sample. Finally we blasted to find other virus genes that contain the most most differentially expressed protein from the extracted cds. 
Along with the python file I have also provided some sample data which is a smaller version of the original donor data. 
With the python script you will also need these packages:
- Biopython
- Kallisto
- Sleuth
- Pandas

When attempting to run the pipeline, please clone the repository. Then please run the python file "PipelineProject.py", it will copy the other files itself and no other action is needed. One db file is missing from the repo as I am unable to upload it to Git due to file size. Instead please download blast+ the run the command makeblastdb -in sequence.fasta -out dbHerpes -title dbHerpes -dbtype nucl this makes the db on your local device using the sequence.fasta file I've provided.  
