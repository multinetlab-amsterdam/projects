## The longitudinal relation between executive functioning and multilayer network topology in glioma patients
### Marike R van Lingen, Lucas C Breedt, Jeroen J G Geurts, Arjan Hillebrand, Martin Klein, Mathilde C M Kouwenhoven, Shanna D Kulik, Jaap C Reijneveld, Cornelis J Stam, Philip C De Witt Hamer, Mona L M Zimmermann, Fernando A N Santos, Linda Douw

Scripts used in doi: 10.1007/s11682-023-00770-w.

<br>

### General Explanations


Scripts: 
	There were two main scripts used for this project:
  
  1. supra_adjacencymatrix_github.m --> to create one supra adjacency matrix outputfile in Matlab including all subjects and their two timepoints. 
  This script follows several steps, for which additional scripts are needed:
    1. fft_filt_BNA.m      --> to filter the MEG signal into six frequencies 
    2. pli_matteo.m        --> to construct the functional network adjacency matrix based on the phase lag index (PLI)
    3. kruskal_algorithm.m --> to construct the minimum spanning tree (MST)
    4. Brain Connectivity Toolbox --> https://www.nitrc.org/frs/?group_id=241 version 2019/03/03
    
  
  2. newest_multilayer_main_code.py --> to calculate multilayer eigenvector centrality in the frontoparietal network, using the supra adjacency matrix output file   	from the supra_adjacencymatrix_github.m script as input file. 
  See the original main multilayer script here: https://github.com/multinetlab-amsterdam/data_analysis/tree/Multilayer/Multilayer
  We set the amount of layers to n = 6 (for each frequency band) and selected the 12 regions from the frontoparietal network as based on the AAL atlas and use the **group_eigenvector_centrality** function:
  
 _function_output(group_eigenvector_centrality, supra_mst,'/mnt/resource/m.vanlingen/m2b/matlab_output_final_okt2021/whole_multilayer_EC', 'EC', list(range(6)))_
   
Statistics:
  The final output from the two abovementioned scripts was one average multilayer eigenvector centrality frontoparietal value per subject per timepoint.
  This was entered in SPSS and statistics were performed as decribed in the method section of the publication. 


</br>

