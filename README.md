# dana

Creation dana project


Notes d'Abigail 11/2/2016
Probleme de merge partiel: 
1. choix du seuil de chevauchement minimal (m nucleotides)
2. si chevauchement(c) inferieur a m (0<c<m), remplacer les c bases du gauche du read R par des X
3. si c == 0, mettre 2 N sur la fin du read F et le debut du read R
4. dans les deux cas, coller les deux bouts

(N et X sont pareils mais comme ca on pourra differencier entre les cas 2 et 3)

Probleme de mismatch dans la partie chevauchement:
Comme le F et le R viennent du meme molecule, on devrait pas avoir des differences, donc on cherche les matches exactes.
On mettra donc pas un 2eme seuil.
On veut jeter la sequence s'il y a des differences dans la partie chevauchee.

Exemple que tu m'as envoye:
ACTATAACAGACGTATTCACAGTAGGATTAAATATCTCAAGTGCGCAAGGCGTAGACGTCGGGCG
ACTATAACAGACGTATTCACAGTAGGATTAAAGATCGAAAGTGGGAAAGGAGTAG

Dans ce cas la, il faut jeter le sequence et faire un note dans le rapport dont tu m'as parle.
