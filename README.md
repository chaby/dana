# dana
======

##1 Merge Process :

Considering 2 DNA fastq sequences:
* the first one a "classic" DNA sequence, read left to right
* the second one a complemented and reversed sequence, read right to left

### 1-1 Convert the complemented and reversed sequence into the "right" way.
### 1-2 apply drastic thresold
  * Cut all the characters at the end of the first sequence with a quality score lesser than the threshold.
  * Cut all the characters at the begining of the second sequence with a quality score lesser than the threshold.

### 1-3 Search the 10 first character of the second sequence into the first sequence
  * Search from the end to the begiging of the first sequence
  * if this 10 characters are found into the first sequence : try to merge the 2 sequences



## 2 The merge problem :


Notes d'Abigail 11/2/2016
Probleme de merge partiel: 
>1. choix du seuil de chevauchement minimal (m nucleotides)
>2. si chevauchement(c) inferieur a m (0<c<m), remplacer les c bases du gauche du read R par des X
>3. si c == 0, mettre 2 N sur la fin du read F et le debut du read R
>4. dans les deux cas, coller les deux bouts

(N et X sont pareils mais comme ca on pourra differencier entre les cas 2 et 3)

Probleme de mismatch dans la partie chevauchement:
Comme le F et le R viennent du meme molecule, on devrait pas avoir des differences, donc on cherche les matches exactes.
On mettra donc pas un 2eme seuil.
On veut jeter la sequence s'il y a des differences dans la partie chevauchee.

Exemple que tu m'as envoye:
ACTATAACAGACGTATTCACAGTAGGATTAAATATCTCAAGTGCGCAAGGCGTAGACGTCGGGCG
ACTATAACAGACGTATTCACAGTAGGATTAAAGATCGAAAGTGGGAAAGGAGTAG

Dans ce cas la, il faut jeter le sequence et faire un note dans le rapport dont tu m'as parle.

