## The relationship between pathological brain activity and functional network connectivity in glioma patients
### Mona LM Zimmermann, Lucas C Breedt, Eduarda GZ Centeno, Jaap C Reijneveld, Fernando AN Santos, Cornelis J Stam, Marike R van Lingen, Menno M Schoonheim, Arjan Hillebrand, Linda Douw

Scripts used in [paper link] 

<br>

### General Explanations


Scripts: 
	There were several different scripts used for this project.

	The scripts used for data preparation are (in the order of use): 

	1. PLI.py: In this script we filter the MEG signal into the different frequency bands and construct the functional network matrix (PLI)  
	2. EC_CC_pipeline: This script is used to threshold and binarize the networks, and then calculate the network metrics (CC and EC) 
	3. masks.py: In this script we define the peritumoral area based om the tumor masks 
	4. master_df.py: This script serves to construct a main dataframe for the main analyses containing the standardized activity and network metrics 
	and the information on which regions belong to which area (peritumoral, homologue, rest of the brain). Based on this dataframe, the area specific dataframes were
	constructed by filtering for the specific area (peritumoral, homologue, rest). This script also
	includes the construction of the dataframes used in the spin test analysis, where raw values are needed. 

	To calculate the offset, we used the scripts found here: https://github.com/multinetlab-amsterdam/data_analysis/tree/master/fooof


	The scripts used for the statistical analyses are: 

	5. statistics_profiles_github.py: This script contains all analyses to characterise activity and network metrics.  
	6. statistics_network_activity_github.py: This script contains all statistical analyses (except for the spin test) to
	relate activity and network metrics.
	7. spin-test-binom.py: This script was used to gain information for the calculation of the binomial 
	confidence interval for the p-values of the spin test results.  


	For the spin test we used the scripts found here: https://github.com/multinetlab-amsterdam/projects/blob/master/LENS_paper_2022/SpinTest.m

For privacy reasons we cannot share all the data used in this project.
</br>
