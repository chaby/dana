#This is a program to assess MiSeq count tables from DANA for DIPLOID LOCI
#this only removes 2 reads and calculates ratios between 1 and 2 to genotype individuals.
#make sure the wd contains all of the count table files for the locus but remove all empty count tables

library(dplyr) #load dplyr 

out_table = NULL #initiate empty data frame

files<-list.files() #create list of all files in the working directory

for (i in files) {
  
  t<-read.table(i,header=FALSE) #read file
  t_sort<-arrange(t,desc(V3)) #sort table by total read counts, largest first
  t_ratio1<-t_sort[2,3]/t_sort[1,3] #calculate ratio btwn first and second largest read counts (btwn 0 and 1)
  #t_ratio2<-t_sort[3,3]/t_sort[1,3] #calculate ratio btwn first and third largest read counts (btwn 0 and 1)
  t_totreads<-sum(t$V3) #calculate the total read number in this well
  y<-abbreviate(i,minlength=6) #truncate file name to get the well name
  
  out_table<-rbind(out_table, data.frame(t_sort[1,1],t_sort[2,1],t_sort[1,3],t_sort[2,3],t_totreads, t_ratio1,i,y)) #create table for output and bind to empty frame
  
}


colnames(out_table)<-c("read1","read2","count1","count2","total","ratio1","file","well") #rename columns for graph

z<-arrange(out_table,desc(total)) #sort output table by decreasing total
#head(z) #print first lines of sorted table


#Then use write.table function to output table as txt file for further use
write.csv(z,"C:/Users/Abigail/Desktop/olon147510_diploid_counts.csv")