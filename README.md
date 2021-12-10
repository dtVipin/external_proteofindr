# proteofindr
 Detect true proteins based on motif search against class  
Copyright (C) 2015 by   
Deepti Vipin, Mishra Lab, CCMB ,India   
All rights reserved   
Released under MIT license (see LICENSE.txt)   


  Requirements
  -----------
-  Perl version 5.20
- FIMO ver >= 4.9.1
- OpalServices



#This script  the FIMO tool for detection of proteins based on motifs.#
#The putative proteins are then filtered to predict true proteins by function.#

Files in inputs
-----------
- 'crp0.meme.xml' Fasta file containing motifs in MEME format
- fastafile of searched proteome data as input (here ex bacterial)

Files in log
-----------
- 'log2.txt' file with UniProt identifiers
- 'log3.txt' contains data downloaded from Uniprot
- 'log4.txt' extracted ID, AA length
- 'log5.txt' contains proteins less than 80AA residues
- 'log6.txt' FIMO result of proteins<80AA

Files in output
-----------
'Fimo.txt' output matching motifs from FIMO
'FINAL.txt' predicted putative proteins
