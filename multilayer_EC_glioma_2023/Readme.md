## [The longitudinal relation between executive functioning and multilayer network topology in glioma patients](https://pubmed.ncbi.nlm.nih.gov/37067658/)
### Marike R van Lingen, Lucas C Breedt, Jeroen J G Geurts, Arjan Hillebrand, Martin Klein, Mathilde C M Kouwenhoven, Shanna D Kulik, Jaap C Reijneveld, Cornelis J Stam, Philip C De Witt Hamer, Mona L M Zimmermann, Fernando A N Santos, Linda Douw

<br>

### General Explanations


There were two main scripts used for this project:
  
  1. **supra_adjacencymatrix_github.m** --> to create one supra adjacency matrix outputfile in Matlab including all subjects and their two timepoints. 
  This script follows several steps, for which additional scripts/toolboxes are needed:
    1. fft_filt_BNA.m      --> to filter the MEG signal into six frequencies 
    2. pli_matteo.m        --> to construct the functional network adjacency matrix based on the phase lag index (PLI)
    3. kruskal_algorithm.m --> to construct the minimum spanning tree (MST)
    4. [Brain Connectivity Toolbox](https://www.nitrc.org/frs/?group_id=241) -->  version 2019/03/03
    
  
  2. **adjusted_multilayer_main_code.py** --> to calculate multilayer eigenvector centrality in the frontoparietal network, using the supra adjacency matrix output file   	from the supra_adjacencymatrix_github.m script as input file. 
  See the original [main multilayer script](https://github.com/multinetlab-amsterdam/data_analysis/tree/Multilayer/Multilayer) that was used as a basis: 
  We 1) set the initial settings, i.e. amount of layers to n = 6 (for each frequency band) and 2) selected the 12 regions from the frontoparietal network as based on the AAL atlas and 3) use the **group_eigenvector_centrality** function:
  
  1) Under headers 'SETTINGS' and 'CREATING LAYER TAGS':
  
    layer_size 	= 78 (AAL has 78 regions, so 78 nodes per layer)
    weighted	= False (because we used the MST matrices that are unweighted and set to either 0 or 1)
    filename	= 'supra_mst_full.mat' (output file from supra_adjacencymatrix_github.m script)
    layer_tags	=['0 = pli delta, 1 = pli theta, 2 = pli alpha1, 3 = pli alpha2, 4 = pli beta, 5 = pli gamma']
    just_tags	=['pli_delta', 'pli_theta', 'pli_alpha1', 'pli_alpha2', 'pli_beta', 'pli_gamma'] 
    plot_tags	=['PLI delta', 'PLI theta', 'PLI alpha1', 'PLI alpha2', 'PLI beta', 'PLI gamma'] 

  2) Under header 'OTHER FUNCTIONS' set sub network (FPN in this case):

    sub_net=[5,6,7,16,17,36,44,45,46,55,56,75]  (number represent FPN regions based on AAL atlas minus 1 since python starts at 0)
  
  3) At the very end of the script, define the function_output:
    
    function_output(group_eigenvector_centrality, supra_mst,'/path/to/whole_multilayer_EC', 'EC', list(range(6)))
   
Statistics:
  The final output from the two abovementioned scripts was one average multilayer eigenvector centrality frontoparietal value per subject per timepoint.
  This was entered in SPSS and statistics were performed as decribed in the method section of the publication. 


</br>

