
## Multimodal multilayer network centrality relates to executive functioning
This repository contains all the files associated with this manuscript by Breedt et al., [available as a preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.28.450180v1.full).


### Main codes and data
All main analyses were performed using two scripts. 

Pre-processed neuroimaging data from multiple modalities (dMRI, fMRI, MEG) were converted into single-layer connectivity matrices and multilayer supra-adjacency matrices in MATLAB using *mumo_script_clean_v1.m*, available in the /projects/mumo/script folder together with a readme.

The supra-adjacency matrices were then transferred to Python, where multilayer eigenvector centrality was computed using *Multilayer_Main_code.py*, [found here](https://github.com/multinetlab-amsterdam/data_analysis/tree/Multilayer/Multilayer).

### Using the main codes and data
To replicate the results of this manuscript, use the supra-adjacency matrix provided in /projects/mumo/data, *mumo_manuscript_supra.mat*, and follow these steps carefully:

**1. Compute nodal multilayer eigenvector centrality using *Multilayer_Main_code.py***
* To ensure functionality, check that the correct version packages are being used.
* Replace **line 66** of *Multilayer_Main_code.py* with 
```python
filename = 'mumo_manuscript_supra.mat'
```
* Run the entire code. This will load the data and initialize the necessary functions.
* Now run the following lines of code to obtain a .csv file containing the nodal multilayer eigenvector centrality values per subject. Rename *filename* to be the desired name of (and path to) the output file and *colname* to be the column name of the column containing the eigenvector centrality values in the output file: 
```python
filename = '/path/to/savefile/filename'
colname = 'nodalMultilayerEc'
data = supra_mst
function_output(group_eigenvector_centrality,data,filename,colname,[0,1,2,3,4,5,6,7])
```
**2. Extract and average multilayer eigenvector centrality of the frontoparietal network using *mumo_script_clean_v1.m***
* Load the deleted regions (found in /projects/mumo/data):
```matlab
load('mumo_manuscript_deleted_regions.mat')
```
* Run lines 464-469 to determine the indices of the areas belonging to the FPN.
* Run lines 572-595 to extract and average multilayer eigenvector centrality of the FPN. Make sure to change the path to the Python-generated .csv file, as well as the save location.

### Statistical analyses
Hierarchical multiple regression models were performed in IBM SPSS Statistics to assess the association between single- and multilayer eigenvector centrality of the FPN and executive functioning. Files and syntax are uploaded, see /projects/mumo/data.