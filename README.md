# PipelineProjectCOMP383

Hello and welcome to a pipeline created by Sid Mehrotra!
For Dr. Wheeler and Daniel: I did track 1
For others:
This pipeline focuses on differential expression between 4 transcriptomes. 
There were two donors, each having samples from 2 days and 6 days post infection with Human herpesvirus 5.
In this pipeline we quantify the transcript levels of different sequences and then use the transcript levels to find differentially expressed genes between the two timepoints (2dpi and 6dpi).
We then blasted to find  what other virus genes contain the most differentially expressed protein from the donors we used.
As you run this pipeline, a log file will be output which details some different things and results from the different steps and tools used in this pipeline. 
Speaking of which! This pipeline uses two big tools among some other libraries. 
Tools: 
- Kallisto: To quantify the transcript levels from the sequences of our donors. 
	- http://pachterlab.github.io/kallisto/ 
	- https://pachterlab.github.io/kallisto-sleuth-workshop-2016/
- Sleuth: To compare the experimental conditions and differential expression after quantifying the transcript levels in TPM
	- https://pachterlab.github.io/sleuth/ 
	- https://pachterlab.github.io/kallisto-sleuth-workshop-2016/
	- Libraries:
		- Sleuth (R)
		- dplyr (R)
- Biopython: Used to pull sequences and parse through them
	- http://biopython.org/DIST/docs/tutorial/Tutorial.html
	- Libraries:
		- Entrez
		- SeqIO
		- Seq
		- SeqRecord
		- NCBIWWW

Libraries:
- Pandas for dataframe manipulation

Now that you should hopefully be sorted in terms of what you need to run this pipeline, how do we run it?

First, clone this repository. There are quite a few files here so I'll explain them. 
- Anything that is with "dbHerpes" was the database files created of all the nucleotide sequences from the subfamily of Betaherpesvirinae as in NCBI. 
- pipelineProject.py is our main python file, this is key
- Since I'm a gentleman, I've also included some sample data, anything labeled "sample"... is that. This data is just shortened versions from our four donors. 
- sleuthScript.R is our R script using the sleuth tool. 

To run this thing, just run the pipelineProject.py file. 
From there it will ask if you want to use sample data, if yes, type "sample", if no, type anything else and it should work (I hope).

