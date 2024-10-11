#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistical analyses for the activity EF project

This script contains all statistical analyses for the project on activity and executive functioning in glioma patients. 

"""
__author__ = "Christina Ulrich & Mona Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "23/05/23"   ### Date it was created
__status__ = "Production" 


####################
# Review History   #
####################

# Reviewed by 


## Script for Statistical Analysis ##

####################
# LIBRARIES       #
####################

# Standard imports ###
import numpy as np 
import pandas as pd 
from scipy import stats
from scipy.stats import shapiro, ttest_rel, wilcoxon
import statsmodels.api as sm
import os
import pyreadstat
import matplotlib.pyplot as plt
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.compat import lzip
import statsmodels.stats.api as sms
import seaborn as sns
# Third party imports ###
# Internal imports ### 



#%%

####################
# DATAFRAME IMPORT #
####################

#ATTENTION>>> to review this code it is easiest to first run the create_df_masks_overlaps.py 
#script ,for the patient data, and then use the created variables for the statistical analyses below.
#For the healthy control data we can use the dataframes previously created by Christina, see below.

#%%
### HCs Hemispheric FPN ###

#HCs (N=25)  
df_HC_Left = pd.read_csv("/data/anw/anw-work/MULTINET/culrich/03_analysis/03_Dataframes/average_HC_Left_df.csv")
df_HC_Right = pd.read_csv("/data/anw/anw-work/MULTINET/culrich/03_analysis/03_Dataframes/average_HC_Right_df.csv")

### HC Info (i.e. info on EF scores and covariates)
HC_info = pd.read_csv('/data/anw/anw-work/MULTINET/culrich/02_participant info/02_HC info/HC_subject_info_after_matching.csv')
#Prepare HC_info df to match other dfs
HC_info.rename(columns = {"Case_ID": "sub"}, inplace = True)
HC_info['sub'] = HC_info['sub'].astype(int).astype(str)
HC_info['sub'] = HC_info['sub'].str.replace('110+', '') #edit sub_ID (remove leading 110) 
HC_info['sub'] = HC_info['sub'].astype(int) #change sub back into int to match other dfs

### Dataframe with all HC data (left and right FPN) ###
#HCs_std = pd.read_csv("/data/anw/anw-work/MULTINET/culrich/03_analysis/02_standardization_for_df/df_activity_standardized_HC.csv") #all regions 

#%%

#%%
### Patient info and EF scores ###
patient_info = pyreadstat.read_sav("/data/anw/anw-work/MULTINET/culrich/02_participant info/01_patient info/Glioomproject_multilayer_BLFU_elekta_n37_forMona_Christina_1.sav")[0]
patient_info['Case_ID'] = 'sub-' + patient_info['Case_ID'].astype(int).astype(str).str.zfill(4)
# Rename colums of EF_score in patient_info file so that both CST and WFT names match"
patient_info.rename(columns = {'flu_dier_corrected_1_zscore' : 'flu_dier_corrected_1_Zscore',
                                'flu_dier_corrected_2_zscore' : 'flu_dier_corrected_2_Zscore' }, inplace = True )


### Lateralization 
li_subs = pd.read_csv("/data/anw/anw-work/MULTINET/m.zimmermann/01_projects/2023_activity_EF/02_analysis/01_subjects/all_subs_list.csv")
#%%

####################
# --- ttest to analyze changes in activity (T1-->T2) --- #
####################

def paired_ttest (df_combined, activity_metrics, area, alpha = 0.05 ):
    '''
    Paramters
    ---------
    df_combined : pd.DataFrame,
        dataframe containing data on activityfor both baseline and FU 
    activity_metrics: str,
        activity that should be investigated (i.e. bbp or offset)
    area: str, 
        determines the area that is being investigated (i.e. peritumoral, ipislateral FPN, contralateral FPN) 
    alpha : float 
        set to 0.05
    '''
    
    #Remove patients from DFs that have progression before FU
    #filter patient_info df based on progression before FU
    patient_info_filtered = patient_info[patient_info['progression_before_FU'] == 2]
    print(patient_info_filtered.shape)
    print(len(patient_info[patient_info['progression_before_FU'] == 1]))
   
    
    df_combined_filt = df_combined[df_combined['sub'].isin(patient_info_filtered['Case_ID'])]
    print(df_combined_filt)
    
       
    # calculate differences between the paired T1 and T2 measurements
    differences = df_combined_filt[f"{activity_metrics}_T2"] - df_combined_filt[f"{activity_metrics}_T1"]
    
    
  
    #testing normality of differnces 
    stat, p = stats.shapiro(differences)
    
    if p > alpha: # null hypothesis (differences has a normal distribution), use paried t-test
        print ("The null hypothesis cannot be rejected --> x has normal distribution")
        print ("Use paired t-test")
        
        #stat_t, p_t = stats.ttest_rel(df_baseline_filtered[activity_metrics], df_FU_filtered[activity_metrics])
        stat_t, p_t = stats.ttest_rel(df_combined_filt[f"{activity_metrics}_T1"], df_combined_filt[f"{activity_metrics}_T2"])
        if p_t < alpha:
                print('Paired t-test: Reject H0 --> there is a significant difference')
        else:
                print ('Paired t-test: Fail to reject H0 --> there is no significant difference')
        print ('T1 vs T2 t = ', str(round(stat_t, 2)), ' p_t = ', str(round(p_t, 4)))
        
        
        
        
    else: # Differences are not normally distributed, use Wilcoxon signed rank test
        print ("The null hypothesis can be rejected --> x is not normally distributed")
        print ("Use Wilcoxon signed-rank test")
        
        #stat_w, p_w = stats.wilcoxon(df_baseline_filtered[activity_metrics], df_FU_filtered[activity_metrics])
        stat_w, p_w = stats.wilcoxon(df_combined_filt[f"{activity_metrics}_T1"],df_combined_filt[f"{activity_metrics}_T2"])
        if p_w < alpha:
            print ('Wilcoxon signed rank test: Reject H0 --> there is a significant difference')
        else:
            print ('Wilcoxon signed rank test: Fail to reject H0 --> no significant difference')
        print ('T1 vs T2 t = ', str(round(stat_w, 2)), ' p = ', str(round(p_w, 4)))
        
       
#%%
#%%    
####################
### Patient paired t-test analysis to investigate change between T1 and T2 ###
####################

#Peritumoral area
paired_ttest(df_peri_combined, 'offset_z', 'peritumoral', alpha = 0.05)
paired_ttest(df_peri_combined,  'BB_welch_z', 'peritumoral', alpha = 0.05)

#%%
#Peri-resection area
paired_ttest(df_cavity_combined, 'offset_z', 'peritumoral', alpha = 0.05)
paired_ttest(df_cavity_combined,  'BB_welch_z', 'peritumoral', alpha = 0.05)

#%%
#Ipsilateral FPN 

result_ipsi_offset = paired_ttest(df_ipsi_combined, 'offset_z', 'ipsilateral', alpha = 0.05)
result_ipsi_bbp = paired_ttest(df_ipsi_combined,'BB_welch_z', 'ipsilateral', alpha = 0.05)
#%%

#Contralateral FPN 
result_contra_offset = paired_ttest(df_contra_combined, 'offset_z', 'contralateral', alpha = 0.05)
result_contra_bbp = paired_ttest(df_contra_combined, 'BB_welch_z', 'contralateral', alpha = 0.05)
    
#%%
####################
### Post-hoc patient paired t-test analysis for the different freq bands ###
####################
#load in frequency specific data

df_peri_baseline_freqs = pd.read_csv('/data/anw/anw-work/MULTINET/m.zimmermann/01_projects/2023_activity_EF/02_analysis/03_dataframes/20240517_std_rel_power_baseline.csv')
df_peri_FU_freqs = pd.read_csv('/data/anw/anw-work/MULTINET/m.zimmermann/01_projects/2023_activity_EF/02_analysis/03_dataframes/20240517_std_rel_power_FU_all.csv')

df_peri_combined_freqs = pd.merge(df_peri_baseline_freqs, df_peri_FU_freqs, on='sub', suffixes=('_T1', '_T2'))


#%%

#paired t-test for the various frequency bands
result_alpha1 = paired_ttest(df_peri_combined_freqs, 'alpha1_z', 'peritumor', alpha = 0.05)

#%%
result_alpha2 = paired_ttest(df_peri_combined_freqs, 'alpha2_z', 'peritumor', alpha = 0.05)

#%%%
result_beta = paired_ttest(df_peri_combined_freqs, 'beta_z', 'peritumor', alpha = 0.05)

#%%
result_delta = paired_ttest(df_peri_combined_freqs, 'delta_z', 'peritumor',alpha = 0.05)

#%%
result_gamma = paired_ttest(df_peri_combined_freqs, 'gamma_z', 'peritumor', alpha = 0.05)

#%%
result_theta = paired_ttest(df_peri_combined_freqs, 'theta_z', 'peritumor',  alpha = 0.05)

#%%
#### Test to test whether activity in the FPN is higher or lower than in HCs ####

def test_activity(df_patients, df_HCs_left, df_HCs_right,df_lateralization, area, activity, time):
    """
    Function to test for a difference in activity between the ipsilateral or contralateral and HC hemispheres. 
    Always matches with the right lateralization (i.e. right ipsilateral matches with right HCs hemisphere).
    Uses the normal unpaired t-test when normality and equality of variances is confirmed. Uses the Welch test when the normality 
    is assumed but variances are unequal. If normality is violated uses the Mann-Whitney U test.  

    Parameters
    ----------
    df_patients : pd.DataFrame,
        dataframe containing all patient data
    df_HCs_left : pd.Dataframe,
        HC dataframe with only left FPN regions
    df_HCs_right : pd.DataFrame,
        HC dataframe with only right FPN regions
    df_lateralization : pd.DataFrame,
        dataframe contiaining info on tumor lateralization (left or right)
    area : str,
        ipsilateral or contralateral: do we compare the ipsilaeral or contralateral hemisphere to HCs?
    activity : str,
        offset_z or BB_welch_z: do we look at offset or bbp?
    time : str,
        baseline or FU, do we look at the baseline measurement or FU?

    Returns
    -------
    pd.DataFrame, 
        results dataframe 

    """
    
   

    df_lat = pd.merge(df_patients, df_lateralization, left_on = "sub", right_on = "Case_ID")
       
    right_tumor = df_lat['lateralization'] == 'right'
    left_tumor = df_lat['lateralization'] == 'left'
    
    #filter for subjects that have a right/left tumor
    df_right = df_lat[right_tumor]
    df_left = df_lat[left_tumor]
    
       
    if area ==  'ipsilateral':
    
        ### Statistical tests ### 
        
        ### Right ipsilateral vs. Right FPN HCs
        results_df_right = ind_ttest(df_right, df_HCs_right,area, activity, 'right_ipsi (i.e. right tumor hemisphere)' )
        
        ### Left ipsilateral vs. Left FPN HCs
        results_df_left = ind_ttest(df_left, df_HCs_left, area, activity, 'left_ipsi (i.e. left tumor hemisphere)')
        
        #results_all = pd.concat([results_df_right,results_df_left])
        
    elif area == 'contralateral':
        
        ### Statistical tests ### 
        
        ### Right ipsilateral vs. Right FPN HCs
        results_df_right= ind_ttest(df_right, df_HCs_left,area, activity, 'right_contra (i.e. left non-tumor hemisphere)' )
        
        ### Left ipsilateral vs. Left FPN HCs
        results_df_left = ind_ttest(df_left, df_HCs_right,area, activity, 'left contra (i.e. right non-tumor hemisphere)')
        
        #results_all = pd.concat([results_df_right,results_df_left])
        
    
    #results_all.to_csv(f'/data/anw/anw-work/MULTINET/m.zimmermann/01_projects/2023_activity_EF/02_analysis/04_results/test_against_HCs_{activity}_{area}_{time}.csv')
   
    #return(results_all)
        
        
#%%

def ind_ttest(df_pat, df_HCs, area, activity, lateralization):
    """
    Function to test for differences with the correct test (after checking the assumptions of the data).
    Used in function test_activity. 

    Parameters
    ----------
    df_pat : pd.DataFrame
        dataframe containing patient data.
    df_HCs : pd.DataFrame, 
        dataframe containing HCs data.
    area : str,
        ipsilateral or contralateral 
    activity : str,
       offset_z or BB_welch_z
    lateralization : str, 
        indicates which area of patients is tested against HCs

    Returns
    -------
    None.

    """
    
    results_df = pd.DataFrame(columns= ['Area', 'activity', 'Lateralization', 'Test', 'T-value', 'pvalue'])
    
    ### Caclculate mean and std of patients and HCs ###
    mean_pat = df_pat[activity].mean()
    std_pat = df_pat[activity].std()
    
    mean_HC = df_HCs[activity].mean()
    std_HC = df_HCs[activity].std()
   
    #testing normality of data
    stat_pat, p_pat = stats.shapiro(df_pat[activity])
    stat_HCs, p_HCs = stats.shapiro(df_HCs[activity])
    
    if (p_pat <0.05) or (p_HCs<0.05): 
        print('Assumption of normality violated: Doing Mann Whitney U test')
        
        results = stats.mannwhitneyu(df_pat[activity], df_HCs[activity])
        print(results)
        
        
        
        # #append results to the results DataFrame
        # #results_df= results_df.append({'Area': area,
        # #                            'activity': activity, 
        #                             'Lateralization': lateralization,
        #                             'Test': 'Mann-Whitney U test, non-normality',
        #                             'U-value':results[0], 
        #                             'pvalue':results[1],
        #                             'mean_pat':mean_pat,
        #                             'std_pat':std_pat,
        #                             'mean_HCs':mean_HC,
        #                             'std_HCs':std_HC},
        #                             ignore_index=True)

    else: 
        print('Assumption of normality not violated: testing for unequal variances now')
       #Testing for equal variances of the two samples (with unequal sample size)
        res = stats.levene(df_pat[activity], df_HCs[activity])
        
        if res.pvalue >0.5: 
            print('Samples have equal variances, doing normal t-test')
            results = stats.ttest_ind(a=df_pat[activity], b=df_HCs[activity], equal_var=True)
            print(results)
            
            # #append results to the results DataFrame
            # results_df= results_df.append({'Area': area,
            #                                 'activity': activity, 
            #                                 'Lateralization': lateralization,
            #                                 'Test': 'Independent t-test, equal_var',
            #                                 'T-value':results[0], 
            #                                 'pvalue':results[1],
            #                                 'mean_pat':mean_pat,
            #                                 'std_pat':std_pat,
            #                                 'mean_HCs':mean_HC,
            #                                 'std_HCs':std_HC},
            #                                 ignore_index=True)
        
        else:
            print('Samples have unequal variances, doing Welchs test')
            results = stats.ttest_ind(a=df_pat[activity], b=df_HCs[activity], equal_var=False)
            print(results)
            # results_df = results_df.append({'Area': area,
            #                                 'activity': activity, 
            #                                 'Lateralization': lateralization,
            #                                 'Test': 'Independent t-test, unequal_var',
            #                                 'T-value':results[0], 
            #                                 'pvalue':results[1],
            #                                 'mean_pat':mean_pat,
            #                                 'std_pat':std_pat,
            #                                 'mean_HCs':mean_HC,
            #                                 'std_HCs':std_HC},
            #                                 ignore_index=True)
        
   
   
  
        

 #%%
### Testing activity against HCs ### 

### Offset_z
#Baseline
test_activity(df_ipsi_averaged, df_HC_Left, df_HC_Right, li_subs,  'ipsilateral',  'offset_z', 'baseline')  
test_activity(df_contra_averaged, df_HC_Left, df_HC_Right, li_subs,  'contralateral',  'offset_z', 'baseline') 

#FU
test_activity(df_ipsi_averaged_FU, df_HC_Left, df_HC_Right, li_subs,  'ipsilateral',  'offset_z', 'FU')  
test_activity(df_contra_averaged_FU, df_HC_Left, df_HC_Right, li_subs,  'contralateral',  'offset_z', 'FU') 

### BBP
#Baseline
test_activity(df_ipsi_averaged, df_HC_Left, df_HC_Right, li_subs,  'ipsilateral',  'BB_welch_z', 'baseline')  
test_activity(df_contra_averaged, df_HC_Left, df_HC_Right, li_subs,  'contralateral',  'BB_welch_z', 'baseline') 
#FU
test_activity(df_ipsi_averaged_FU, df_HC_Left, df_HC_Right, li_subs,  'ipsilateral',  'BB_welch_z', 'FU')  
test_activity(df_contra_averaged_FU, df_HC_Left, df_HC_Right, li_subs,  'contralateral',  'BB_welch_z', 'FU')    




#%%
####################
# --- Linear regressions --- #  --> Multiple Linear Regression
####################

####################
### 1) to test relationship Activity vs EF (T1) ###
####################
def run_multiple_linear_regression(df_activity, df_ef, EF_Score, area, output_dir, covariates=[]):
    """
    Parameters:
    
    df_activity: pandas.DataFrame,
        the dataframe containing the activity_metrics (bbp and offset) 
    df_ef: pandas.DataFrame,
        the dataframe containing the patient characteristics (EF scores and covariates)
    EF_Score: str,
        the column name in df representing the dependent variable (i.e. EF Z scores: either CST or WFT)
    area: str, 
          determines the area that is being investigated (i.e. preitumoral, ipsilateral FPN, contralateral FPN)
    output_dir: str, 
          path to directory where results should be stored
    covariates: list of str,
          a list of column names in df representing the covariates (default: [])
    """
  
    li = []     
         
    for activity_metrics in ["BB_welch_z", "offset_z"]: 
              
        # Merge the two dataframes based on sub_ID
        merge_df = pd.merge(df_activity, df_ef, left_on='sub', right_on='Case_ID')
        print(merge_df)    

             
        #In case EF_Score = CST is analzyed --> remove sub-9038 and sub-0087 from DF
        if EF_Score == "cstc_corrected_1_Zscore":
            merge_df.drop(merge_df[(merge_df["sub"] == "sub-9038")|(merge_df["sub"] == "sub-0087")].index, inplace = True)
            merge_df.reset_index(inplace = True)
            
            #Check to see if correct participant was excluded
            print(merge_df["sub"])   
            print(merge_df["Case_ID"])
            
        else:
            #Check to see if correct participant are still included for WFT
            print(merge_df["sub"])   
            print(merge_df["Case_ID"])
            
        
        # Define the dependent variable (outcome) and independent variables
        X_cols = [activity_metrics] + covariates
        print(X_cols)
        
        X = merge_df[X_cols]
        y = merge_df[EF_Score]
        
        #Check if X inludes corect columns 
        print(X)
        #Check if Y inludes corect columns 
        print(y)
        
        
        # Add a constant term to the independent variables
        X = sm.add_constant(X)

        # Create a multiple linear regression model
        model = sm.OLS(y, X)

        # Fit the model to the data
        results = model.fit()

        # Print the regression coefficients and other results
        print(f"Results for {activity_metrics} vs {EF_Score} in {area}:")
        print (results.summary())
     
       

      
       ### --- Testing Assumption of Multiple Linear Regression --- ###
        print ('--- Testing Assumptions of Multiple Linear Regression ---')
        
        #Extract the independent variables (IV) only
        X_IV = X.drop(columns=['const'])
        print(X_IV.head())
        
        # 1) Linear Relationship (between Indepentent and target variables)
        print('1) Testing for Linear Relationship: See Graphs')
        for IVs in X_IV.columns:
              plt.scatter(X[IVs], y)
              plt.xlabel(IVs)
              plt.ylabel(EF_Score)
              plt.title(f"Scatter Plot of {IVs} vs {EF_Score}")
              plt.show()
            
        # 2) No Multicollinearity
        print ('2) Testing for Multicollinearity')
          
        #Caculate the correlation matrix
        corr_matrix = X_IV.corr()
        print(corr_matrix)
        
        # Check for variables with high correlation coefficients
        high_corr = set()
        for i in range(len(corr_matrix.columns)):
            for j in range(i):
                if abs(corr_matrix.iloc[i, j])> 0.7:
                    colname = corr_matrix.columns[i]
                    high_corr.add(colname)
        print("Variables with high correlation coefficients:", high_corr)
        
        #Calculate VIF scores for each independent variable (IVs)
        
        vif_scores = pd.DataFrame()
        vif_scores["feature"] = X_IV.columns
        vif_scores["VIF"] = [variance_inflation_factor(X_IV.values, i) for i in range(X_IV.shape[1])]
        print(vif_scores)
        

        #3) Homoscedasticity - constant variance
        print ("3) Testing Homoscedasticity: See Graph")
        residuals = results.resid
        fitted_vals=results.predict(X)
        plt.scatter(fitted_vals, residuals)
        plt.xlabel('Fitted Values')
        plt.ylabel('Residuals')
        plt.title("Scatter Plot to test Homoscedasticity")
        plt.show()
        
        name = ['Lagranage multiplier statistic', 'p-value', 'f-value', 'f p-value']
        test = sms.het_breuschpagan(results.resid, results.model.exog)
        print(lzip(name, test))
        
        
        #4) No Autocorrelation of errors 
        print ("4) Testing Autocorrelation of errors: See Graph")
        plt.plot(residuals.index, residuals)
        plt.title("Plot to test Autocorrelation of errors")
        plt.show()
        
        #5) Residual Normality 
        print("5) Testing Residual Normality: See Graph")
        #sns.distplot(residuals).set(title ="Residual plot: Testing Normality")
        
        print('residual mean:')
        print(np.mean(residuals))
        
        stat, p = stats.shapiro(results.resid)
        print('Shapiro test:')
        print(stat, p)
        
        #6) Checking for outliers
        print("6) Checking for Residual relation with independent variables/ outliers: See Graph")
        fig, axs = plt.subplots(ncols=X_IV.shape[1], figsize=(15, 5))
        for i, col in enumerate (X_IV.columns):
            axs[i].scatter(X_IV[col], results.resid)
            axs[i].set_xlabel(col)
            axs[i].set_ylabel("Residuals")
            plt.title("Checking for outliers")
        plt.show()
        
        #or
        residuals =results.resid
        sm.qqplot(residuals, line='s')
        plt.title('Residuals Q-Q')
        plt.show()
        
        #Calculate Cook's distance values 
        influence = results.get_influence()
        cooks = influence.cooks_distance[0]
        
        #Plot the Cook's distance values against the obersvation numbers
        plt.plot(cooks, 'o')
        plt.xlabel('Observation number')
        plt.ylabel("Cook's distance")
        plt.title("Observation number vs Cook's distance")
        plt.show()
        
        #Identify influential observations 
        # n= len(y)
        # threshold = 4/n
        # influential_observations = np.where(cooks >= threshold)[0]
        # print('Influential observations:', influential_observations)
        
        
        print("--- Testing Assumptions of Multiple Linear Regression Completed --- ")
        
     
        # get the summary table of the results
        table1 = results.summary2().tables[0]
        table2 = results.summary2().tables[1]
        table2.reset_index(inplace=True)
      
        summary_df = pd.concat([table1, table2], axis=1)
        summary_df = summary_df.reset_index()
      #Append results to the dataframe
        li.append(summary_df)
                
    results_df = pd.concat(li, axis=0) #indent seems not to matter ? (22-05-24)

    
    # Save the dataframe to CSV
    results_df.to_csv(f'{output_dir}linear_regression_baseline_{area}_{EF_Score}.csv')
    #return(results_df)
#%%

###############################
### Patient linear regression_baseline analysis ### 
###############################
    output_dir = '/data/anw/anw-work/MULTINET/m.zimmermann/01_projects/2023_activity_EF/02_analysis/04_results/'
         
#Peritumoral
    #CST
    run_multiple_linear_regression(df_peri_averaged, patient_info, 'cstc_corrected_1_Zscore', 'peritumoral', output_dir, covariates=["Dummy_location_frontal_or_not", 'Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted'])#,"epilepsy_dich.1"])
    #2023-06-28 MZ: epilepsy_dich taken out due to multicollinearity (VIF >0.5)
#%%
    #WFT
    run_multiple_linear_regression(df_peri_averaged, patient_info, 'flu_dier_corrected_1_Zscore', 'peritumoral', output_dir, covariates=["Dummy_location_frontal_or_not", 'Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted'])
   

#%%                                    
#Ipsilateral FPN
    #CST
    run_multiple_linear_regression(df_ipsi_averaged, patient_info, 'cstc_corrected_1_Zscore', 'ipsilateral', output_dir, covariates =["Dummy_location_frontal_or_not", 'Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted'])#, "epilepsy_dich.1"])
    #2023-06-28 MZ: epilepsy_dich taken out due to multicollinearity (VIF >0.5)
#%%
    
    #WFT
    run_multiple_linear_regression(df_ipsi_averaged, patient_info, 'flu_dier_corrected_1_Zscore', 'ipsilateral', output_dir, covariates=["Dummy_location_frontal_or_not", 'Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted'])
#%%                                              
#Contralateral FPN
    #CST#
    run_multiple_linear_regression(df_contra_averaged, patient_info, 'cstc_corrected_1_Zscore', 'contralateral', output_dir, covariates=["Dummy_location_frontal_or_not", 'Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted'])#, "epilepsy_dich.1"])
    #2023-06-28 MZ: epilepsy_dich taken out due to multicollinearity (VIF >0.5)

#%%

    #WFT
    run_multiple_linear_regression(df_contra_averaged, patient_info, 'flu_dier_corrected_1_Zscore', 'contralateral', output_dir, covariates=["Dummy_location_frontal_or_not", 'Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted'])


    
#%%
####################   
### 2) to test relationship Activity vs EF (longitudinally) ###
####################

   
def run_delta_regression(df_activity, df_ef, EF_Score, area, output_dir, covariates=[]):
    """
  Parameters:
    df_activity: pandas.DataFrame,
        the dataframe containing the activity_metrics (bbp and offset) for both timepoints (baseline and FU)
    df_ef: pandas.DataFrame,
        the dataframe containing the patient characteristics (EF scores and covariates)
    EF_Score: str,
        the column name in df representing the dependent variable (i.e. EF Z scores: either CST or WFT)
    area: str, 
          determines the area that is being investigated (i.e. preitumoral, ipsilateral FPN, contralateral FPN)
    output_dir: str, 
          path to directory where results should be stored
    covariates: list of str,
          a list of column names in df representing the covariates (default: [])
    """
 

    li = []              
    for activity_metrics in ["BB_welch_z", "offset_z"]: 
        
        # Merge the two dataframes based on sub_ID
        merge_df_org = pd.merge(df_activity, df_ef, left_on='sub', right_on='Case_ID')
        merge_df = pd.merge(df_activity, df_ef, left_on='sub', right_on='Case_ID')

        #In case EF_Score = CST is analzyed --> remove sub-9038 and sub-0087 from DF
        if EF_Score == "cstc_corrected":
            merge_df.drop(merge_df[(merge_df["sub"] == "sub-9038")|(merge_df["sub"] == "sub-0087")].index, inplace = True)
            merge_df.reset_index(inplace = True)
            #Check to see if correct participant was excluded
            print(merge_df["sub"])   
            print(merge_df["Case_ID"])
            
            
        #In case EF_Score = WFT is analyzed --> remove sub-9031 from DF
        else:
            #Check to see if correct participant are still included for WFT
            merge_df.drop(merge_df[(merge_df["sub"] == "sub-9031")].index, inplace = True)
            merge_df.reset_index(inplace = True)
            #Check to see if correct participant was excluded
            print(merge_df["sub"])   
            print(merge_df["Case_ID"])
        
        
        #Remove patients from DF that have progression before FU
        prog_filt = merge_df["progression_before_FU"] == 1
        merge_df = merge_df[~prog_filt]
        
        merge_df.reset_index(inplace = True)
   
        #Check to see if correct participants and correct columns are included in the final version    
        print(merge_df["sub"])  
        print(merge_df.columns)
    
        
        # Calculate delta values
        merge_df[f'delta_{activity_metrics}'] = merge_df[activity_metrics + '_T2'] - merge_df[activity_metrics + '_T1']
      
        merge_df[f'delta_{EF_Score}']= merge_df[EF_Score + "_2_Zscore"] - merge_df[EF_Score + "_1_Zscore"]
       
       
    
        # Define the dependent variable (outcome) and independent variables
        if len(covariates) > 0:
            X = pd.concat([merge_df[f'delta_{activity_metrics}'], merge_df[covariates]], axis=1)
        else:
            X = merge_df[f'delta_{activity_metrics}']
        print(X)  
        
        Y = merge_df[f'delta_{EF_Score}']
        print(Y)

        #Add a constant term to the independent variables
        X = sm.add_constant(X)
        
        # Create a multiple linear regression model 
        model = sm.OLS(Y, X)
        #Fit the model to the data
        results = model.fit()

        # Print the regression coefficients and other results
        print(f"Results for delta {activity_metrics} vs delta {EF_Score} in {area}:")
        print (results.summary())
        
      
    
       ### --- Testing Assumption of Multiple Linear Regression --- ###
        print ('--- Testing Assumptions of Multiple Linear Regression ---')
        
        #Extract the independent variables (IV) only
        X_IV = X.drop(columns=['const'])
        print(X_IV.head())
        
        # 1) Linear Relationship (between Indepentent and target variables)
        print('1) Testing for Linear Relationship: See Graphs')
        for IVs in X_IV.columns:
              plt.scatter(X[IVs], Y)
              plt.xlabel(IVs)
              plt.ylabel(EF_Score)
              plt.title(f"Scatter Plot of {IVs} vs {EF_Score}")
              plt.show()
            
        # 2) No Multicollinearity
        print ('2) Testing for Multicollinearity')
          
        #Caculate the correlation matrix
        corr_matrix = X_IV.corr()
        print(corr_matrix)
        
        # Check for variables with high correlation coefficients
        high_corr = set()
        for i in range(len(corr_matrix.columns)):
            for j in range(i):
                if abs(corr_matrix.iloc[i, j])> 0.7:
                    colname = corr_matrix.columns[i]
                    high_corr.add(colname)
        print("Variables with high correlation coefficients:", high_corr)
        
        #Calculate VIF scores for each independent variable (IVs)
        
        vif_scores = pd.DataFrame()
        vif_scores["feature"] = X_IV.columns
        vif_scores["VIF"] = [variance_inflation_factor(X_IV.values, i) for i in range(X_IV.shape[1])]
        print(vif_scores)
        

        #3) Homoscedasticity - constant variance
        print ("3) Testing Homoscedasticity: See Graph")
        residuals = results.resid
        fitted_vals=results.predict(X)
        plt.scatter(fitted_vals, residuals)
        plt.xlabel('Fitted Values')
        plt.ylabel('Residuals')
        plt.plot(Y, [0]*len(Y))
        plt.title("Scatter Plot to test Homoscedasticity")
        plt.show()
        
        name = ['Lagranage multiplier statistic', 'p-value', 'f-value', 'f p-value']
        test = sms.het_breuschpagan(results.resid, results.model.exog)
        print(lzip(name, test))
        
        #4) No Autocorrelation of errors 
        print ("4) Testing Autocorrelation of errors: See Grahp")
        plt.plot(residuals.index, residuals)
        plt.title("Plot to test Autocorrelation of errors")
        plt.show()
        
        #5) Residual Normality 
        print("5) Testing Residual Normality: See Graph")
        sns.distplot(residuals).set(title ="Residual plot: Testing Normality")
        print('residual mean:')
        #print(np.mean(residuals))
        
        stat, p = stats.shapiro(results.resid)
        print('Shapiro test:')
        print(stat, p)
        
        #6) Checking for outliers
        print("6) Checking for Residual relation with independent variables/ outliers: See Graph")        
        fig, axs = plt.subplots(ncols=X_IV.shape[1], figsize=(15, 5))
        for i, col in enumerate (X_IV.columns):
            axs[i].scatter(X_IV[col], results.resid)
            axs[i].set_xlabel(col)
            axs[i].set_ylabel("Residuals")
            plt.title("Checking for outliers")
        plt.show()
        
        #or
        residuals =results.resid
        sm.qqplot(residuals, line='s')
        plt.title('Residuals Q-Q')
        plt.show()
        
        print("--- Testing Assumptions of Multiple Linear Regression Completed --- ")
          
        # get the summary table of the results
        table1 = results.summary2().tables[0]
        table2 = results.summary2().tables[1]
        table2.reset_index(inplace=True)
      
        summary_df = pd.concat([table1, table2], axis=1)
        summary_df = summary_df.reset_index()
      #Append results to the dataframe
        li.append(summary_df)
                
    results_delta_df = pd.concat(li, axis=0)
        

    # Save the dataframe to CSV
    results_delta_df.to_csv(f'{output_dir}/linear_regression_delta_scores_{area}_{EF_Score}.csv')
    #return(results_delta_df)
    

#%%
###############################
### Patient linear regression_delta_scores analysis ### 
###############################
    output_dir = "/data/anw/anw-work/MULTINET/m.zimmermann/01_projects/2023_activity_EF/02_analysis/04_results"
         
#Peritumoral
    #CST
    run_delta_regression(df_peri_combined, patient_info, 'cstc_corrected', 'peritumoral', output_dir, covariates=['Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted',"Dummy_CT_during_FU"])
    #2023-06-28 MZ: Covariate "Interval_surgery_NPA" taken out after multicollinearity
#%%

    #WFT
    run_delta_regression(df_peri_combined, patient_info, 'flu_dier_corrected', 'peritumoral', output_dir, covariates=['Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted','Dummy_RT', 'Dummy_RTandXT'])
    
#%%                                
#Ipsilateral FPN
    #CST
    run_delta_regression(df_ipsi_combined, patient_info, 'cstc_corrected', 'ipsilateral', output_dir, covariates=['Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted',"Dummy_CT_during_FU", "Interval_surgery_NPA"])
  
#%%    
  #WFT
    run_delta_regression(df_ipsi_combined, patient_info, 'flu_dier_corrected', 'ipsilateral', output_dir, covariates=['Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted','Dummy_RT', 'Dummy_RTandXT'])
#%%                                                
#Contralateral FPN
    #CST
    run_delta_regression(df_contra_combined, patient_info, 'cstc_corrected', 'contralateral', output_dir, covariates=['Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted',"Dummy_CT_during_FU", "Interval_surgery_NPA"])
        #Interval_Surgery_MEG taken out due to multicolinearity# NOT when redid analysis (28-06-23)!
#%%        
        #WFT
    run_delta_regression(df_contra_combined, patient_info, 'flu_dier_corrected', 'contralateral', output_dir, covariates=['Dummy_IDH_WT', 'Dummy_IDHmut_noncodeleted','Dummy_RT', 'Dummy_RTandXT'])

#%%          
####################
### 3) to test relationship Activity vs EF in HCs ###
####################

def run_HC_linear_regression(df_activity, df_ef, EF_Score, area, output_dir, covariates=[]):
    """
    Parameters:
    df_activity: pandas.DataFrame,
        the dataframe containing the activity_metrics (bbp and offset) 
    df_ef: pandas.DataFrame,
        the dataframe containing the patient characteristics (EF scores and covariates)
    EF_Score: str,
        the column name in df representing the dependent variable (i.e. EF Z scores: either CST or WFT)
    area: str, 
          determines the area that is being investigated (i.e. left or right FPN)
    output_dir: str, 
          path to directory where results should be stored
    covariates: list of str,
          a list of column names in df_ef representing the covariates (default: [])
    """
  
    li = []
              
    for activity_metrics in ["BB_welch_z", "offset_z"]: 
              
        # Merge the two dataframes based on sub_ID
        merged_df = pd.merge(df_activity, df_ef, on='sub')
        print(merged_df.head)
        print(merged_df[[activity_metrics, EF_Score]])
        
        # Define the dependent variable (outcome) and independent variables
        X_cols = [activity_metrics] + covariates
        X = merged_df[X_cols]
        y = merged_df[EF_Score]

        #Check if X inludes corect columns 
        print(X)
        #Check if Y inludes corect columns 
        print(y)
 
        # Add a constant term to the independent variables
        X = sm.add_constant(X)

        # Create a multiple linear regression model
        model = sm.OLS(y, X)

        # Fit the model to the data
        results = model.fit()

        # Print the regression coefficients and other results
        print(f"Results for {activity_metrics} vs {EF_Score} in {area}:")
        print (results.summary())
   
       ### --- Creating Graphs --- ###
        sns.lmplot(x=X.columns[1], y=EF_Score, data=merged_df, scatter_kws={'alpha':0.5})
        plt.title(f'for HCs: {activity_metrics} vs. {EF_Score}_in {area}')
        plt.axhline(y= -1.5, linestyle ='--', color ='grey')
      #  plt.savefig(f'/data/anw/anw-work/MULTINET/culrich/03_analysis/05_Graphs/04_HC_Linear_Regression/Plot of {activity_metrics} vs {EF_Score} in HCs in {area}".png', bbox_inches = "tight")
        plt.show()
        
 
    ### --- Testing Assumption of Multiple Linear Regression --- ###
        print ('--- Testing Assumptions of inear Regression ---')
         
        #Extract the independent variables (IV) only
        X_IV = X.drop(columns=['const'])
        print(X_IV.head())
         
        # 1) Linear Relationship (between Indepentent and target variables)
        print('1) Testing for Linear Relationship: See Graphs')
        for IVs in X_IV.columns:
             plt.scatter(X[IVs], y)
             plt.xlabel(IVs)
             plt.ylabel(EF_Score)
             plt.title(f"Scatter Plot of {IVs} vs {EF_Score}")
             plt.show()
            
        #2) Homoscedasticity - constant variance
        print ("2) Testing Homoscedasticity: See Graph")
        residuals = results.resid
        fitted_vals=results.predict(X)
        plt.scatter(fitted_vals, residuals)
        plt.xlabel('Fitted Values')
        plt.ylabel('Residuals')
        plt.plot(y, [0]*len(y))
        plt.title("Scatter Plot to test Homoscedasticity")
        plt.show()
        
        name = ['Lagranage multiplier statistic', 'p-value', 'f-value', 'f p-value']
        test = sms.het_breuschpagan(results.resid, results.model.exog)
        print(lzip(name, test))
        
        #3) Independence of observations
        #create plot of the residuals against the independent variable
        print('3) Testing for Independence: See Graphs')
        for IVs in X_IV.columns:
            plt.scatter(X[IVs], residuals)
            plt.xlabel(IVs)
            plt.ylabel('residuals')
            plt.title(f"Scatter Plot of {IVs} vs residuals")
            plt.show()
            
        #4) No Autocorrelation of errors 
        print ("4) Testing Autocorrelation of errors: See Grahp")
        plt.plot(residuals.index, residuals)
        plt.title("Plot to test Autocorrelation of errors")
        plt.show()
        
        #5) Residual Normality 
        print("5) Testing Residual Normality: See Graph")
        sns.distplot(residuals).set(title ="Residual plot: Testing Normality")
        
        print('residual mean:')
        print(np.mean(residuals))
        
        stat, p = stats.shapiro(results.resid)
        print('Shapiro test:')
        print(stat, p)
        
        #6) Checking for outliers
        print("6) Checking for outliers: See Graph")
    
        residuals =results.resid
        sm.qqplot(residuals, line='s')
        plt.title('Residuals Q-Q')
        plt.show()
        
        print("--- Testing Assumptions of Linear Regression Completed --- ")
        
      # get the summary table of the results
        table1 = results.summary2().tables[0]
        table2 = results.summary2().tables[1]
        table2.reset_index(inplace=True)
        
        summary_df = pd.concat([table1, table2], axis=1)
        summary_df = summary_df.reset_index()
        li.append(summary_df)
              
    results_HC_df = pd.concat(li, axis=0)

    # Save the dataframe to CSV
    results_HC_df.to_csv(output_dir + f'/linear_regression_HC_{area}_{EF_Score}.csv')
    #return(results_HC_df)
   
#%%
###############################
### HC linear regression_baseline analysis ### 
###############################
output_dir = '/data/anw/anw-work/MULTINET/culrich/03_analysis/04_Output_Statistics/04_HC_Linear_Regression'
         
#Left Hemisphere FPN
    #CST
run_HC_linear_regression(df_HC_Left, HC_info, 'CST_Zscore_shifting', 'left_FPN', output_dir, covariates=[])
    #WFT
run_HC_linear_regression(df_HC_Left, HC_info, 'AnimalFluency_Zscore', 'left_FPN', output_dir, covariates=[])
                      
 #Rigt Hemisphere FPN
    #CST
run_HC_linear_regression(df_HC_Right, HC_info, 'CST_Zscore_shifting', 'right_FPN', output_dir, covariates=[])
    #WFT
run_HC_linear_regression(df_HC_Right, HC_info, 'AnimalFluency_Zscore', 'right_FPN', output_dir, covariates=[])
 
