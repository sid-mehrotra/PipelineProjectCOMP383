import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio.Blast import NCBIWWW
sampleOrReal = input("If you would like to run sample data, please type 'sample': ")
os.makedirs("PipelineProject_Sid_Mehrotra") #making a new directory where everything will be stored
os.system("cp sleuthScript.R PipelineProject_Sid_Mehrotra")
os.system("cp dbHerpes.ndb PipelineProject_Sid_Mehrotra")
os.system("cp dbHerpes.nhr PipelineProject_Sid_Mehrotra")
os.system("cp dbHerpes.nin PipelineProject_Sid_Mehrotra")
os.system("cp dbHerpes.not PipelineProject_Sid_Mehrotra")
os.system("cp dbHerpes.nsq PipelineProject_Sid_Mehrotra")
os.system("cp dbHerpes.ntf PipelineProject_Sid_Mehrotra")
os.system("cp dbHerpes.nto PipelineProject_Sid_Mehrotra")
if sampleOrReal == "sample":
	os.system("cp sampleSRR5660030_1.fastq PipelineProject_Sid_Mehrotra")
	os.system("cp sampleSRR5660030_2.fastq PipelineProject_Sid_Mehrotra")
	os.system("cp sampleSRR5660033_1.fastq PipelineProject_Sid_Mehrotra")
	os.system("cp sampleSRR5660033_2.fastq PipelineProject_Sid_Mehrotra")
	os.system("cp sampleSRR5660044_1.fastq PipelineProject_Sid_Mehrotra")
	os.system("cp sampleSRR5660044_2.fastq PipelineProject_Sid_Mehrotra")
	os.system("cp sampleSRR5660045_1.fastq PipelineProject_Sid_Mehrotra")
	os.system("cp sampleSRR5660045_2.fastq PipelineProject_Sid_Mehrotra")
	os.chdir("PipelineProject_Sid_Mehrotra") #go into the directory we just made
else:
	os.chdir("PipelineProject_Sid_Mehrotra") #go into the directory we just made
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030'") #obtain the 4 transcriptomes
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033'")
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044'")
	os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045'")

	os.system("fastq-dump -I --split-files SRR5660030") #uncompress the data you downloaded
	os.system("fastq-dump -I --split-files SRR5660033")
	os.system("fastq-dump -I --split-files SRR5660044")
	os.system("fastq-dump -I --split-files SRR5660045")


Entrez.email = "smehrotra1@luc.edu"
handle = Entrez.efetch(db = "nucleotide", id = " NC_006273.2", rettype = "gb", retmode = "text") #found from biopython tutorial book thing, grabbing the genome in genbank form
record = SeqIO.read(handle, format = "genbank") #reads it nicely and separates it

#now to get the CDS, we can use the features and add it to one file https://www.biostars.org/p/9542737/#9547242

numCDS = 0 #keeping track of num of coding sequences
cds = [] #array keeping track of all the coding sequwnces
mySeqRecords = list()
outfile = open("fastaFile.fasta", "a")
for f in record.features: #for every feature
	if f.type == "CDS": #if it is a CDS
		#print(f)
		cds.append(f.extract(record.seq)) #append the sequences, got this from https://widdowquinn.github.io/2018-03-06-ibioic/01-introduction/03-parsing.html
		numCDS += 1 #add 1 to our counter
		mySeq = f.extract(record.seq)
		mySeqRecord = SeqRecord(mySeq)
		#print(f.qualifiers["protein_id"])  #you can access all the stuff because f is technically a dictionary
					#https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
		mySeqRecord.id = f.qualifiers["protein_id"]
		#print(mySeqRecord)
		#print(type(mySeqRecords))
		mySeqRecords.append(mySeqRecord)
		stringToWrite = ">" + str(mySeqRecord.id)
		outfile.write(stringToWrite + "\n")
		outfile.write(str(mySeqRecord.seq) + "\n")

		#print(len(mySeqRecord.seq))
outfile.close()

logFile = open("PipelineProject.log", "a")
logFile.write("The HCMV genome (NC_006273.2) has " + str(numCDS) + " CDS. \n" + "\n")
logFile.close()

os.system("time kallisto index -i index.idx fastaFile.fasta")
if sampleOrReal == "sample":
	os.makedirs("results/SRR5660030") #making a new directory where everything will be stored
	os.system("time kallisto quant -i index.idx -o results/SRR5660030 -b30 -t 2 sampleSRR5660030_1.fastq sampleSRR5660030_2.fastq")
	os.makedirs("results/SRR5660033")
	os.system("time kallisto quant -i index.idx -o results/SRR5660033 -b30 -t 2 sampleSRR5660033_1.fastq sampleSRR5660033_2.fastq")
	os.makedirs("results/SRR5660044")
	os.system("time kallisto quant -i index.idx -o results/SRR5660044 -b30 -t 2 sampleSRR5660044_1.fastq sampleSRR5660044_2.fastq")
	os.makedirs("results/SRR5660045")
	os.system("time kallisto quant -i index.idx -o results/SRR5660045 -b30 -t 2 sampleSRR5660045_1.fastq sampleSRR5660045_2.fastq")
else:
	os.makedirs("results/SRR5660030") #making a new directory where everything will be stored
	os.system("time kallisto quant -i index.idx -o results/SRR5660030 -b30 -t 2 SRR5660030_1.fastq SRR5660030_2.fastq")
	os.makedirs("results/SRR5660033")
	os.system("time kallisto quant -i index.idx -o results/SRR5660033 -b30 -t 2 SRR5660033_1.fastq SRR5660033_2.fastq")
	os.makedirs("results/SRR5660044")
	os.system("time kallisto quant -i index.idx -o results/SRR5660044 -b30 -t 2 SRR5660044_1.fastq SRR5660044_2.fastq")
	os.makedirs("results/SRR5660045")
	os.system("time kallisto quant -i index.idx -o results/SRR5660045 -b30 -t 2 SRR5660045_1.fastq SRR5660045_2.fastq")

# Passing the TSV file to
# read_csv() function
# with tab separator
# This function will
# read data from file
SRR5660030KallistoResults = pd.read_csv('results/SRR5660030/abundance.tsv', sep='\t')
SRR5660033KallistoResults = pd.read_csv('results/SRR5660033/abundance.tsv', sep='\t')
SRR5660044KallistoResults = pd.read_csv('results/SRR5660044/abundance.tsv', sep='\t')
SRR5660045KallistoResults = pd.read_csv('results/SRR5660045/abundance.tsv', sep='\t')
# printing data
#print(SRR5660030KallistoResults)
with open('PipelineProject.log', 'a') as f:
	f.write("sample" + "\t" + "condition" + "\t" + "min_tpm" + "\t" + "med_tpm" + "\t" + "mean_tpm" + "\t" + "max_tpm" + "\n")
	f.write("Donor 1" + "\t" + "2dpi" + "\t" + str(SRR5660030KallistoResults["tpm"].min()) + "\t" + str(SRR5660030KallistoResults["tpm"].median()) + "\t" + str(SRR5660030KallistoResults["tpm"].mean()) + "\t" + str(SRR5660030KallistoResults["tpm"].max()) + "\n")
	f.write("Donor 1" + "\t" + "6dpi" + "\t" + str(SRR5660033KallistoResults["tpm"].min()) + "\t" + str(SRR5660033KallistoResults["tpm"].median()) + "\t" + str(SRR5660033KallistoResults["tpm"].mean()) + "\t" + str(SRR5660033KallistoResults["tpm"].max()) + "\n")

	f.write("Donor 2" + "\t" + "2dpi" + "\t" + str(SRR5660044KallistoResults["tpm"].min()) + "\t" + str(SRR5660044KallistoResults["tpm"].median()) + "\t" + str(SRR5660044KallistoResults["tpm"].mean()) + "\t" + str(SRR5660044KallistoResults["tpm"].max()) + "\n")
	f.write("Donor 2" + "\t" + "6dpi" + "\t" + str(SRR5660045KallistoResults["tpm"].min()) + "\t" + str(SRR5660045KallistoResults["tpm"].median()) + "\t" + str(SRR5660045KallistoResults["tpm"].mean()) + "\t" + str(SRR5660045KallistoResults["tpm"].max()) + "\n" + "\n")


f.close()

with open("sleuth_table.txt", "a") as f:
	f.write("sample" + "\t" +  "condition" + "\t" + "path" + "\n")
	f.write("SRR5660030" + "\t" + "2dpi" + "\t" + "results/SRR5660030" + "\n")
	f.write("SRR5660033" + "\t" + "6dpi" + "\t" + "results/SRR5660033" + "\n")
	f.write("SRR5660044" + "\t" + "2dpi" + "\t" + "results/SRR5660044" + "\n")
	f.write("SRR5660045" + "\t" + "6dpi" + "\t" + "results/SRR5660045" + "\n")
f.close()

os.system("Rscript sleuthScript.R")

sleuthResults = pd.read_csv("fdr_results.txt", sep = " ")
with open("PipelineProject.log", "a") as f:
	f.write("target_id" + "\t" + "test_stat" + "\t" + "pval" + "\t" + "qval" + "\n")
	for i in range(0, len(sleuthResults)):
		f.write(str(sleuthResults.iloc[i,0]) + "\t" + str(sleuthResults.iloc[i,3]) + "\t" + str(sleuthResults.iloc[i, 1]) + "\t" + str(sleuthResults.iloc[i, 2]) + "\n")
f.close()

topHit = str(sleuthResults.iloc[0,0])
print(topHit)
handle = Entrez.efetch(db = "protein", id = "YP_081530.1", rettype = "fasta", retmode = "text")
print(handle)
record = SeqIO.read(handle, format = "fasta")
print(record)

SeqIO.write(record, "TopHitProtein.fasta", "fasta")

input_file = "TopHitProtein.fasta"
output_file = "blastResults.csv"
blast_command = "tblastn -query " + input_file + " -db dbHerpes -out " + output_file + " -outfmt '10 sacc pident length qstart qend sstart send bitscore evaluue stitle'"
os.system(blast_command)

blastResults = pd.read_csv("blastResults.csv", on_bad_lines='skip')
print(blastResults)
with open("PipelineProject.log", "a") as f:
	f.write("sacc" + "\t" + "pident" + "\t" + "length" + "\t" + "qstart" + "\t" + "qend" + "\t" + "sstart" + "\t" + "send" + "\t" + "bitscore" + "\t" + "evalue" + "\t" + "stitle"  +"\n")
	for i in range(0, 10):
		f.write(str(blastResults.iloc[i,0]) + "\t" + str(blastResults.iloc[i,1]) + "\t" + str(blastResults.iloc[i,2]) + "\t" + str(blastResults.iloc[i,3]) + "\t" + str(blastResults.iloc[i,4]) + "\t" + str(blastResults.iloc[i,5]) + "\t" + str(blastResults.iloc[i,6]) + "\t" + str(blastResults.iloc[i,7]) + "\t" + str(blastResults.iloc[i,8]) + "\t" + str(blastResults.iloc[i,9]) + "\n")
f.close()
