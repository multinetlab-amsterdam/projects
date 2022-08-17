The scripts posted here were used in [Modeling of individual neurophysiological brain connectivity](https://www.biorxiv.org/content/10.1101/2022.03.02.482608v1).

For privacy reasons, we cannot share the raw data alongisde the scripts. Therefore, our aim here was to post a general version of the computations, allowing readers to understand how we computed the different parts of the analysis.

Step 1: Running the model
- run_model.m: Run the Jansen-Rit model with the individual structural connectivity as input. 

Step 2: Correlations between simulations and empirical data
- corr_sim_mod.m: Calculate the correlations between simulated and empirical data for all functional connectivity metrics. Subsequently find the optimal coupling value which gives the highest correlation between simulated and empirical data. 
- average_SC_emp_all.m: Calculate the correlations between simulated data with the average structural connectivity as input and empirical data for the AEC long and AEC over epochs. 

Step 3: Higher coupling ranges
- check_higher_coupling_ranges.m: Check the coupling value per subject if the coupling value is the maximal value of the range this subject has to be ran again with the higher coupling range. Subsequently calculate whether the correlation between simulation and empirical data is higher compared to the coupling from the lower range.
- test_between_maximum_coupling_values.m: Test whether there is a difference between the maximal correlations for all FC metrics. 
- test_FC_measures_max_correlations.m: Test whether there is a difference between the maximal coupling values per FC metric.

Step 4: Non-matched versus matched correlations
- random_coupling.m: Comparison between the maximum correlations of subjects own data with the maximum correlations of simulated data with the empirical data of other subjects. 

Step 5: Creating raincloud figures
-	Different raincloud figures were created to compare the maximal correlations, optimal coupling values and to compare the maximal correlations when using the individual structural connectivity with the maximal correlations when using the average structural connectivity. 
