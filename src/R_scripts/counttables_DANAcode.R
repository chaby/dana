#This is a program to assess MiSeq count tables based on DANA code
#To be run in a working directory containing only count tables 
#NB: will break if there are empty count tables in the folder

library(dplyr) #load dplyr 

out_table = NULL #initiate empty data frame

files<-list.files() #create list of all files in the working directory
x<-length(files) #calculate the number of files in the working directory

for (i in files) {
  
  t<-read.table(i,header=FALSE) #read file
  t_sort<-arrange(t,desc(V3)) #sort table by total read counts, largest first
  t_ratio<-t_sort[2,3]/t_sort[1,3] #calculate ratio btwn first and second largest read counts (btwn 0 and 1)
  t_totreads<-sum(t$V3) #calculate the total read number in this well
  y<-abbreviate(i,minlength=6) #truncate file name to get the well name
  
  out_table<-rbind(out_table, data.frame(t_sort[1,1],t_sort[1,3],t_ratio,t_totreads,i,y)) #create table for output and bind to empty frame
}


colnames(out_table)<-c("read","count","ratio","total","file","well") #rename columns for graph

plot(out_table$well,out_table$ratio, ylim=c(0,1),xlab=("well"),ylab=("Ratio")) #graph of the ratio variable for each well

z<-arrange(out_table,desc(ratio)) #sort output table by decreasing ratio
head(z) #print first lines of sorted table

#Then use write.csv function to output table as csv file for further use