###This code will identify thresholds for trimming sequences at a specified thresholds.
###Need to be in a working directory with files that have the well name, the read name (2 columns), and the length.

lengthfiles<-list.files() #create list of all files in the working directory

outputs<-data.frame() #initiate data frame

x0<-0 #set start of bp sequence
breakpoints<-seq(x0, 250, 1) #create vector from 1 to 250

for (i in lengthfiles) {
  
  t<-read.table(i,header=FALSE) #read csv
  size<- sort(breakpoints, decreasing=TRUE) #sort the vector of possible bp in descending order
  total_reads<-length(t$V4) #find total number of reads for the marker
  threshold<-0.9*total_reads #calculate a threshold; 90% of reads will be this long or longer
  value <-c() #initiate blank vector
  marker <-c() #initiate blank vector
  
  
  for (bp in size) {
    if (sum(t$V4>=bp) > threshold) {  #find the number of reads that are this long or longer; compare to threshold
      value <-c(value, bp) #make vector of values that cross the threshold
      
    } 
    
  }
  
  markeroutput<-cbind(i,value[1],max(t$V4)) #make dataframe to connect marker name with the highest value that passes threshold and the max length
  outputs<-rbind(outputs,markeroutput) #bind dataframe for each marker
}

colnames(outputs)<-c("marker","cut_size","max") #add column names
print(outputs) #print dataframe for all markers