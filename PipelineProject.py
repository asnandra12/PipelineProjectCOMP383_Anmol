import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
import pandas as pd
svr = input("Please choose which data to run with the program, for sample please enter 'sample' else hit enter: ")# Here we choose whether to use the sample data provided or obtain the transcriptomes via wget  
os.makedirs("PipelineProject_Anmol_Singh") #making a new directory where everything will be stored
os.system("copy sleuthScript.R PipelineProject_Anmol_Singh") #copy rscript and db files into new directory
os.system("copy dbHerpes.ndb PipelineProject_Anmol_Singh")#note this is run on windows so for it to work the command is copy
os.system("copy dbHerpes.nhr PipelineProject_Anmol_Singh")
os.system("copy dbHerpes.nin PipelineProject_Anmol_Singh")
os.system("copy dbHerpes.not PipelineProject_Anmol_Singh")
os.system("copy dbHerpes.nsq PipelineProject_Anmol_Singh")
os.system("copy dbHerpes.ntf PipelineProject_Anmol_Singh")
os.system("copy dbHerpes.nto PipelineProject_Anmol_Singh")
if svr == "sample":
	os.system("copy sampleSRR5660030_1.fastq PipelineProject_Anmol_Singh")#copy sample data
	os.system("copy sampleSRR5660030_2.fastq PipelineProject_Anmol_Singh")
	os.system("copy sampleSRR5660033_1.fastq PipelineProject_Anmol_Singh")
	os.system("copy sampleSRR5660033_2.fastq PipelineProject_Anmol_Singh")
	os.system("copy sampleSRR5660044_1.fastq PipelineProject_Anmol_Singh")
	os.system("copy sampleSRR5660044_2.fastq PipelineProject_Anmol_Singh")
	os.system("copy sampleSRR5660045_1.fastq PipelineProject_Anmol_Singh")
	os.system("copy sampleSRR5660045_2.fastq PipelineProject_Anmol_Singh")
	os.chdir("PipelineProject_Anmol_Singh") #move into the directory we just made
else: #if sample was not chosen the wget will grab the required files
	os.chdir("PipelineProject_Anmol_Singh")#move into new directory
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030'") #obtain the 4 transcriptomes
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033'")
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044'")
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045'")

	os.system("fastq-dump -I --split-files SRR5660030") #Unpack the downloaded data into the directory
	os.system("fastq-dump -I --split-files SRR5660033")
	os.system("fastq-dump -I --split-files SRR5660044")
	os.system("fastq-dump -I --split-files SRR5660045")


Entrez.email = "asingh16@luc.edu"
handle = Entrez.efetch(db = "nucleotide", id = " NC_006273.2", rettype = "gb", retmode = "text") #pulled from biopython review notes
record = SeqIO.read(handle, format = "genbank") #read and seperate genome

numCDS = 0 #tracking the # of cds
seqnum = 1 #Making sure all seq in the index file are unique
cds = [] #an array for tracking all the sequences
mySeqRecords = list()
outfile = open("fastaFile.fasta", "a")#In order to pull the cds we must look at the features, I used https://www.biostars.org/p/83058/ and https://www.biostars.org/p/9542737/ as a ref to start the code
for f in record.features: #for every feature
	if f.type == "CDS": #if it is a CDS
		cds.append(f.extract(record.seq)) #append the sequences
		numCDS += 1 #updating the counter var to keep track of the cds
		mySeq = f.extract(record.seq)
		mySeqRecord = SeqRecord(mySeq)
		mySeqRecord.id = f.qualifiers["protein_id"]
		mySeqRecords.append(mySeqRecord)#write out the sequences to the fasta file
		stringToWrite = ">" + str(mySeqRecord.id) + str(seqnum)#important to add seqnum here as kallisto won't run with duplicate headers
		outfile.write(stringToWrite + "\n")
		outfile.write(str(mySeqRecord.seq) + "\n")
		seqnum += 1 #updating sequence# to make sure all seq in file are unique and runs fine for kallisto
outfile.close()

logFile = open("PipelineProject.log", "a")
logFile.write("The HCMV genome (NC_006273.2) has " + str(numCDS) + " CDS. \n" + "\n")# write out the results to the log file
logFile.close()

os.system(" kallisto index -i index.idx fastaFile.fasta") #make the index for kallisto quantification from fasta file
if svr == "sample": #using sample files
	os.makedirs("results/SRR5660030") #making a directory for results
#quantification
	os.system(" kallisto quant -i index.idx -o results/SRR5660030 -b30 -t 2 sampleSRR5660030_1.fastq sampleSRR5660030_2.fastq")#note- in the notes I saw that the time command was also attached before kallisto, I'm running this in windows and it doesn't run the same plus I don't understand the use so I removed it
	os.makedirs("results/SRR5660033")
	os.system(" kallisto quant -i index.idx -o results/SRR5660033 -b30 -t 2 sampleSRR5660033_1.fastq sampleSRR5660033_2.fastq")
	os.makedirs("results/SRR5660044")
	os.system(" kallisto quant -i index.idx -o results/SRR5660044 -b30 -t 2 sampleSRR5660044_1.fastq sampleSRR5660044_2.fastq")
	os.makedirs("results/SRR5660045")
	os.system(" kallisto quant -i index.idx -o results/SRR5660045 -b30 -t 2 sampleSRR5660045_1.fastq sampleSRR5660045_2.fastq")
else: #using original files
	os.makedirs("results/SRR5660030") #making a directory for results
	#quantification
	os.system(" kallisto quant -i index.idx -o results/SRR5660030 -b30 -t 2 SRR5660030_1.fastq SRR5660030_2.fastq")
	os.makedirs("results/SRR5660033")
	os.system(" kallisto quant -i index.idx -o results/SRR5660033 -b30 -t 2 SRR5660033_1.fastq SRR5660033_2.fastq")
	os.makedirs("results/SRR5660044")
	os.system(" kallisto quant -i index.idx -o results/SRR5660044 -b30 -t 2 SRR5660044_1.fastq SRR5660044_2.fastq")
	os.makedirs("results/SRR5660045")
	os.system(" kallisto quant -i index.idx -o results/SRR5660045 -b30 -t 2 SRR5660045_1.fastq SRR5660045_2.fastq")

SRR5660030KallistoResults = pd.read_csv('results/SRR5660030/abundance.tsv', sep='\t') #Here we use pandas read_csv to read the tsv results from kallisto
SRR5660033KallistoResults = pd.read_csv('results/SRR5660033/abundance.tsv', sep='\t')
SRR5660044KallistoResults = pd.read_csv('results/SRR5660044/abundance.tsv', sep='\t')
SRR5660045KallistoResults = pd.read_csv('results/SRR5660045/abundance.tsv', sep='\t')

with open('PipelineProject.log', 'a') as f:
	f.write("sample" + "\t" + "condition" + "\t" + "min_tpm" + "\t" + "med_tpm" + "\t" + "mean_tpm" + "\t" + "max_tpm" + "\n")# write our final min, median, mean and max of the results to the log file along with a header
	f.write("Donor 1" + "\t" + "2dpi" + "\t" + str(SRR5660030KallistoResults["tpm"].min()) + "\t" + str(SRR5660030KallistoResults["tpm"].median()) + "\t" + str(SRR5660030KallistoResults["tpm"].mean()) + "\t" + str(SRR5660030KallistoResults["tpm"].max()) + "\n")
	f.write("Donor 1" + "\t" + "6dpi" + "\t" + str(SRR5660033KallistoResults["tpm"].min()) + "\t" + str(SRR5660033KallistoResults["tpm"].median()) + "\t" + str(SRR5660033KallistoResults["tpm"].mean()) + "\t" + str(SRR5660033KallistoResults["tpm"].max()) + "\n")

	f.write("Donor 2" + "\t" + "2dpi" + "\t" + str(SRR5660044KallistoResults["tpm"].min()) + "\t" + str(SRR5660044KallistoResults["tpm"].median()) + "\t" + str(SRR5660044KallistoResults["tpm"].mean()) + "\t" + str(SRR5660044KallistoResults["tpm"].max()) + "\n")
	f.write("Donor 2" + "\t" + "6dpi" + "\t" + str(SRR5660045KallistoResults["tpm"].min()) + "\t" + str(SRR5660045KallistoResults["tpm"].median()) + "\t" + str(SRR5660045KallistoResults["tpm"].mean()) + "\t" + str(SRR5660045KallistoResults["tpm"].max()) + "\n" + "\n")


f.close()

with open("sleuth_table.txt", "a") as f: #write out the table required to run sleuth
	f.write("sample" + "\t" +  "condition" + "\t" + "path" + "\n")#template from notes
	f.write("SRR5660030" + "\t" + "2dpi" + "\t" + "results/SRR5660030" + "\n")
	f.write("SRR5660033" + "\t" + "6dpi" + "\t" + "results/SRR5660033" + "\n")
	f.write("SRR5660044" + "\t" + "2dpi" + "\t" + "results/SRR5660044" + "\n")
	f.write("SRR5660045" + "\t" + "6dpi" + "\t" + "results/SRR5660045" + "\n")
f.close()

os.system("Rscript sleuthScript.R") #Running the slueth script pulled directly from the notes
sleuthResults = pd.read_csv("fdr_results.txt", sep = " ") #read in and wite out sleuth results to log file
with open("PipelineProject.log", "a") as f:
	f.write("target_id" + "\t" + "test_stat" + "\t" + "pval" + "\t" + "qval" + "\n")
	for i in range(0, len(sleuthResults)):
		f.write(str(sleuthResults.iloc[i,0]) + "\t" + str(sleuthResults.iloc[i,3]) + "\t" + str(sleuthResults.iloc[i, 1]) + "\t" + str(sleuthResults.iloc[i, 2]) + "\n")
f.close()

topHit = str(sleuthResults.iloc[0,0]) #finding the top hit sequence 
print(topHit)
handle = Entrez.efetch(db = "protein", id = "YP_081530.1", rettype = "fasta", retmode = "text") #YP_081530.1 is the tophit id obtained from the previous segment
print(handle)
record = SeqIO.read(handle, format = "fasta") 
print(record)

SeqIO.write(record, "TopHitProtein.fasta", "fasta") #writing the fasta info of the TopHit Protein to a file
query_seqfile = "TopHitProtein.fasta" #blast tophit vs Betaherpesvirinae subfamily
output_file = "blastResults.csv"
blast_command = "tblastn -query " + query_seqfile + " -db dbHerpes -out " + output_file + " -outfmt '10 sacc pident length qstart qend sstart send bitscore evaluue stitle'"#blast command from notes with specific parameters added
os.system(blast_command)#note- we use tblastn instead of blastn above due to the fact that the database we are using has many nucleotide sequences inside it

blastResults = pd.read_csv("blastResults.csv", on_bad_lines='skip')#read results using pandas
print(blastResults)
with open("PipelineProject.log", "a") as f: #write out blast results
	f.write("sacc" + "\t" + "pident" + "\t" + "length" + "\t" + "qstart" + "\t" + "qend" + "\t" + "sstart" + "\t" + "send" + "\t" + "bitscore" + "\t" + "evalue" + "\t" + "stitle"  +"\n")
	for i in range(0, 10):#write out the top 10 hits
		f.write(str(blastResults.iloc[i,0]) + "\t" + str(blastResults.iloc[i,1]) + "\t" + str(blastResults.iloc[i,2]) + "\t" + str(blastResults.iloc[i,3]) + "\t" + str(blastResults.iloc[i,4]) + "\t" + str(blastResults.iloc[i,5]) + "\t" + str(blastResults.iloc[i,6]) + "\t" + str(blastResults.iloc[i,7]) + "\t" + str(blastResults.iloc[i,8]) + "\t" + str(blastResults.iloc[i,9]) + "\n")
f.close()
