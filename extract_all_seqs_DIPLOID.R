#This is a program that takes a list of read names and well names and creates a fasta file of the reads we want for a single locus.
#use for DIPLOID LOCI!!!
#Revamped for use after DANA
#input: text file with list of reads and well names associated (2 reads per well!!!)
#for now, header on text file needs to read "read1", "read2", "ratio1", "well"
#also need working directory to contain all fasta files for a locus 
#goodreads table is based on the count_table_DANA_diploid code

library (seqinr) #load seqinr

goodreads<-read.table("C:/Users/Abigail/Desktop/olon_147510_goodreads.txt",header=TRUE) #read in list of names we need; to change each time the file changes
files<-list.files() #make a list of all files in the working directory
seq_out_table_hom = NULL #initialize empty data frame
seq_out_table_het = NULL

hom<-goodreads[goodreads$ratio1<0.3,] #subset all homozygotes

het<-goodreads[goodreads$ratio1>0.6,] #subset all heterozygotes (these two steps will toss some individuals)

for (i in hom$well) {
  namelist<-grep(i,files, value=TRUE) #find file that matches the well (works b/c file names start w/ well names)
  namefasta<-read.fasta(namelist,as.string=TRUE,set.attributes=FALSE) #read the appropriate fasta
  y<-cbind(namefasta[names(namefasta) %in% hom$read1],i) #extract the read in the fasta which matches the read name in the text file, bind with well name
  z<-cbind(namefasta[names(namefasta) %in% hom$read1],i) #extract the first read for the well AGAIN
  
  
  seq_out_table_hom<-rbind(seq_out_table_hom,y,z ) #create output table of HOMOZYGOUS sequences and well names 
  
}


for (i in het$well) {
  namelist<-grep(i,files, value=TRUE) #find file that matches the well (works b/c file names start w/ well names)
  namefasta<-read.fasta(namelist,as.string=TRUE,set.attributes=FALSE) #read the appropriate fasta
  y<-cbind(namefasta[names(namefasta) %in% het$read1],i) #extract the read in the fasta which matches the read name in the text file, bind with well name
  z<-cbind(namefasta[names(namefasta) %in% het$read2],i) #extract the second read for the well
  
  
  seq_out_table_het<-rbind(seq_out_table_het,y,z ) #create output table of HETEROZTGOUS sequences and well names 
  
}

seq_out_table_all<-rbind(seq_out_table_hom,seq_out_table_het) #make table with all individuals

colnames(seq_out_table_all)<-c("read","well") #rename columns on output table; prob not needed

write.fasta(seq_out_table_all[,1],names=seq_out_table_all[,2],"C:/Users/Abigail/Desktop/olon_147510_allout.fasta") 
#command to save to fasta with sequences and well names; change file name for each locus