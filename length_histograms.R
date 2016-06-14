###This code will create histograms to show the distribution of read lengths.
###Need to be in a working directory with files that have the well name, the read name (2 columns), and the length.

lengthfiles<-list.files() #create list of all files in the working directory

x0<-0 #set start of breakpoint sequence
breakpoints<-seq(x0, 250, 1) #create 1 bp bins for histogram

for (i in lengthfiles) {
  
  t<-read.table(i,header=FALSE) #read file
  hist(t$V4,main=i,breaks=breakpoints) #create histogram
  
}

