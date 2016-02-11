# dana

1 Merge Process :
-----------------

Considering 2 DNA fastq sequences:
* the first one a "classic" DNA sequence, read left to right
* the second one a complemented and reversed sequence, read right to left

1-1 Convert the complemented and reversed sequence into the "right" way
1-2 apply drastic thresold
  * Cut all the characters at the end of the first sequence with a quality score lesser than the threshold.
  * Cut all the characters at the begining of the second sequence with a quality score lesser than the threshold.
1-3 Search the 10 first character of the second sequence into the first sequence
  * if this 10 characters are found into the first sequence : try to merge the 2 sequences



2 The merge problem :
-------------------

when we found 

Creation dana project
