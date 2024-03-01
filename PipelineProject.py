import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
import pandas as pd
svr = input("Please choose which data to run with the program, for sample please enter 'sample' else hit enter: ")# Here we choose whether to use the sample data provided or obtain the transcriptomes via wget  
os.makedirs("PipelineProject_Anmol_Singh") #making a new directory where everything will be stored
os.system("copy sleuthScript.R PipelineProject_Anmol_Singh") #copy rscript and db files into new directory
os.system("copy dbHerpes.nhr PipelineProject_Anmol_Singh")#note this is run on windows so for it to work the command is copy
os.system("copy dbHerpes.nin PipelineProject_Anmol_Singh")
os.system("copy dbHerpes.nsq PipelineProject_Anmol_Singh")
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
	os.system("wget 'https://www.ncbi.nlm.nih.gov/sra/SRX2896360'") #obtain the 4 transcriptomes
	os.system("wget 'https://www.ncbi.nlm.nih.gov/sra/SRX2896363'")
	os.system("wget 'https://www.ncbi.nlm.nih.gov/sra/SRX2896374'")
	os.system("wget 'https://www.ncbi.nlm.nih.gov/sra/SRX2896375'")
	os.system("fastq-dump -I --split-files SRR5660030") #Unpack the downloaded data into the directory
	os.system("fastq-dump -I --split-files SRR5660033")
	os.system("fastq-dump -I --split-files SRR5660044")
	os.system("fastq-dump -I --split-files SRR5660045")


Entrez.email = "asingh16@luc.edu"
handle = Entrez.efetch(db = "nucleotide", id = " NC_006273.2", rettype = "gb", retmode = "text") #pulled from biopython review notes, pulling from genebank and recieving text
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

logFile = open("PipelineProject.log", "a")#open the file in append so all written in data will be added to the end
logFile.write("The HCMV genome (NC_006273.2) has " + str(numCDS) + " CDS. \n" + "\n")# write out the results to the log file
logFile.close()

os.system(" kallisto index -i index.idx fastaFile.fasta") #make the index for kallisto quantification from fasta file
if svr == "sample": #using sample files
	os.makedirs("results/SRR5660030") #making a directory for results
	os.system(" kallisto quant -i index.idx -o results/SRR5660030 -b30 -t 2 sampleSRR5660030_1.fastq sampleSRR5660030_2.fastq")#note- in the notes I saw that the time command was also attached before kallisto, I'm running this in windows and it doesn't run the same plus I don't understand the use so I removed it
	os.makedirs("results/SRR5660033")
	os.system(" kallisto quant -i index.idx -o results/SRR5660033 -b30 -t 2 sampleSRR5660033_1.fastq sampleSRR5660033_2.fastq")#quantification
	os.makedirs("results/SRR5660044")
	os.system(" kallisto quant -i index.idx -o results/SRR5660044 -b30 -t 2 sampleSRR5660044_1.fastq sampleSRR5660044_2.fastq")
	os.makedirs("results/SRR5660045")
	os.system(" kallisto quant -i index.idx -o results/SRR5660045 -b30 -t 2 sampleSRR5660045_1.fastq sampleSRR5660045_2.fastq")
else: #using original files
	os.makedirs("results/SRR5660030") #making a directory for results
	os.system(" kallisto quant -i index.idx -o results/SRR5660030 -b30 -t 2 SRR5660030_1.fastq SRR5660030_2.fastq")#quantification
	os.makedirs("results/SRR5660033")
	os.system(" kallisto quant -i index.idx -o results/SRR5660033 -b30 -t 2 SRR5660033_1.fastq SRR5660033_2.fastq")
	os.makedirs("results/SRR5660044")
	os.system(" kallisto quant -i index.idx -o results/SRR5660044 -b30 -t 2 SRR5660044_1.fastq SRR5660044_2.fastq")
	os.makedirs("results/SRR5660045")
	os.system(" kallisto quant -i index.idx -o results/SRR5660045 -b30 -t 2 SRR5660045_1.fastq SRR5660045_2.fastq")

KR_30 = pd.read_csv('results/SRR5660030/abundance.tsv', sep='\t') #Here we use pandas read_csv to read the tsv results from kallisto
KR_33 = pd.read_csv('results/SRR5660033/abundance.tsv', sep='\t')
KR_44 = pd.read_csv('results/SRR5660044/abundance.tsv', sep='\t')
KR_45 = pd.read_csv('results/SRR5660045/abundance.tsv', sep='\t')

KR_30_min = str(KR_30["tpm"].min())# write out final min, median, mean and max of the results
KR_30_median = str(KR_30["tpm"].median())
KR_30_mean = str(KR_30["tpm"].mean())
KR_30_max = str(KR_30["tpm"].max())

KR_33_min = str(KR_33["tpm"].min())# write out final min, median, mean and max of the results
KR_33_median = str(KR_33["tpm"].median())
KR_33_mean = str(KR_33["tpm"].mean())
KR_33_max = str(KR_33["tpm"].max())

KR_44_min = str(KR_44["tpm"].min())# write out final min, median, mean and max of the results
KR_44_median = str(KR_44["tpm"].median())
KR_44_mean = str(KR_44["tpm"].mean())
KR_44_max = str(KR_44["tpm"].max())

KR_45_min = str(KR_45["tpm"].min())# write out final min, median, mean and max of the results
KR_45_median = str(KR_45["tpm"].median())
KR_45_mean = str(KR_45["tpm"].mean())
KR_45_max = str(KR_45["tpm"].max())

with open('PipelineProject.log', 'a') as f:
	f.write("sample" + "\t" + "condition" + "\t" + "min_tpm" + "\t" + "med_tpm" + "\t" + "mean_tpm" + "\t" + "max_tpm" + "\n")# write out final min, median, mean and max of the results to the log file along with a header
	f.write("Donor 1" + "\t" + "2dpi" + "\t" + KR_30_min + "\t" + KR_30_median + "\t" + KR_30_mean + "\t" + KR_30_max + "\n")
	f.write("Donor 1" + "\t" + "6dpi" + "\t" + KR_33_min + "\t" + KR_33_median + "\t" + KR_33_mean + "\t" + KR_33_max + "\n")
	f.write("Donor 2" + "\t" + "2dpi" + "\t" + KR_44_min + "\t" + KR_44_median + "\t" + KR_44_mean + "\t" + KR_44_max + "\n")
	f.write("Donor 2" + "\t" + "6dpi" + "\t" + KR_45_min + "\t" + KR_45_median + "\t" + KR_45_mean + "\t" + KR_45_max + "\n" + "\n")
f.close()

with open("sleuth_table.txt", "a") as f: #write out the table required to run sleuth
	f.write("sample" + "\t" +  "condition" + "\t" + "path" + "\n")#template from notes
	f.write("SRR5660030" + "\t" + "2dpi" + "\t" + "results/SRR5660030" + "\n")
	f.write("SRR5660033" + "\t" + "6dpi" + "\t" + "results/SRR5660033" + "\n")
	f.write("SRR5660044" + "\t" + "2dpi" + "\t" + "results/SRR5660044" + "\n")
	f.write("SRR5660045" + "\t" + "6dpi" + "\t" + "results/SRR5660045" + "\n")
f.close()

os.system("Rscript sleuthScript.R") #Running the slueth script pulled directly from the notes
sleuthResults = pd.read_csv("fdr_results.txt", sep = " ") #read in and write out sleuth results to log file
with open("PipelineProject.log", "a") as f:
	f.write("target_id" + "\t" + "test_stat" + "\t" + "pval" + "\t" + "qval" + "\n")#write out headers
	for i in range(0, len(sleuthResults)):
		f.write(str(sleuthResults.iloc[i,0]) + "\t" + str(sleuthResults.iloc[i,3]) + "\t" + str(sleuthResults.iloc[i, 1]) + "\t" + str(sleuthResults.iloc[i, 2]) + "\n")
f.close()

mdCDS = str(sleuthResults.iloc[0,0]) #finding the most differentially expressed CDC sequence 
handle = Entrez.efetch(db = "protein", id = "YP_081530.1", rettype = "fasta", retmode = "text") #YP_081530.1 is the mdCDS id obtained from the previous segment
record = SeqIO.read(handle, format = "fasta") 
SeqIO.write(record, "MDCDS.fasta", "fasta") #writing the fasta info of the mdCDS Protein to a file
query_seqfile = "MDCDS.fasta" #blast mdCDS vs Betaherpesvirinae subfamily
output_file = "blastResults.csv"
blast_command = "tblastn -query " + query_seqfile + " -db dbHerpes -out " + output_file + " -outfmt '10 sacc pident length qstart qend sstart send bitscore evaluue stitle'"#blast command from notes with specific parameters added
os.system(blast_command)#note- we use tblastn instead of blastn above due to the fact that we are using a protien fasta to blast a nucleotide db
#Using a regualar blastn would not work as it is a nucleotide to nucleotide blast
blastResults = pd.read_csv("blastResults.csv", on_bad_lines='skip')#read results using pandas
with open("PipelineProject.log", "a") as f: #write out blast results
	f.write("sacc" + "\t" + "pident" + "\t" + "length" + "\t" + "qstart" + "\t" + "qend" + "\t" + "sstart" + "\t" + "send" + "\t" + "bitscore" + "\t" + "evalue" + "\t" + "stitle"  +"\n")
	x = 0
	while x < 10 :#write out the top 10 hits
		f.write(str(blastResults.iloc[x,0]) + "\t" + str(blastResults.iloc[x,1]) + "\t" + str(blastResults.iloc[x,2]) + "\t" + str(blastResults.iloc[x,3]) + "\t" + str(blastResults.iloc[x,4]) + "\t" + str(blastResults.iloc[x,5]) + "\t" + str(blastResults.iloc[x,6]) + "\t" + str(blastResults.iloc[x,7]) + "\t" + str(blastResults.iloc[x,8]) + "\t" + str(blastResults.iloc[x,9]) + "\n")
		x += 1
f.close()
