## The longitudinal relation between executive functioning and multilayer network topology in glioma patients
### Marike R van Lingen, Lucas C Breedt, Jeroen J G Geurts, Arjan Hillebrand, Martin Klein, Mathilde C M Kouwenhoven, Shanna D Kulik, Jaap C Reijneveld, Cornelis J Stam, Philip C De Witt Hamer, Mona L M Zimmermann, Fernando A N Santos, Linda Douw

Scripts used in doi: 10.1007/s11682-023-00770-w.

<br>

### General Explanations


Scripts: 
	There were two main scripts used for this project:
  
  1. supra_adjacencymatrix_github.m --> to create one supra adjacency matrix outputfile in Matlab including all subjects and the two timepoints 
  This script follows several steps, for which additional scripts are needed:
    1. fft_filt_BNA.m      --> to filter the MEG signal into six frequencies 
    2. pli_matteo.m        --> to construct the functional network adjacency matrix based on the phase lag index (PLI)
    3. kruskal_algorithm.m --> to construct the minimum spanning tree (MST)
    
    4. adjacency matrix/check matrices script --> described these? Or not used in the end? 
  
  2. newest_multilayer_main_code.py --> to calculate multilayer eigenvector centrality in the frontoparietal network, using the supra adjacency matrix as input file. 
  See the main multilayer script here: https://github.com/multinetlab-amsterdam/data_analysis/tree/Multilayer/Multilayer
  
  
  
Statistics:
  The final output from the two abovementioned scripts was one average eigenvector centrality frontoparietal value per subject per timepoint.
  This was entered in SPSS and statistics were performed as decribed in the method section of the publication. 


</br>

