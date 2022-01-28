* Encoding: UTF-8.
* Analysis accompanying the LENS paper

*** get patient chars

USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.
    
FREQUENCIES VARIABLES=sex histology side graad IDH_1p19q KPS_split hvkbsch
  /ORDER=ANALYSIS.

FREQUENCIES VARIABLES=age tumor_volume_cm3 
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.


USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1 AND progressie = 1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES= PFS  
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.

USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1 AND died = 1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES= OS  
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.


* per subgroup

SORT CASES  BY IDH_1p19q.
SPLIT FILE LAYERED BY IDH_1p19q.

FREQUENCIES VARIABLES=sex histology side graad IDH_1p19q KPS_split hvkbsch
  /ORDER=ANALYSIS.

FREQUENCIES VARIABLES=age tumor_volume_cm3 
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.


USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1 AND progressie = 1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES= PFS  
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.

USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1 AND died = 1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES= OS  
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.


FILTER OFF.
USE ALL.
EXECUTE.

SPLIT FILE OFF.


*** diffs with HCs

T-TEST GROUPS=group(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=age
  /CRITERIA=CI(.95).

CROSSTABS
  /TABLES=group BY sex
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT
  /COUNT ROUND CELLUSE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.


    
********* SURVIVAL ANALYSES with LENS

* IDHmut, codeleted

USE ALL.
COMPUTE filter_$=(IDH_1p19q =0 AND Incl_tumor_probability =1).
VARIABLE LABELS filter_$ 'IDH_1p19q =0 AND Incl_tumor_probability =1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.


COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER bbp_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER slope_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER offset_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* Cox regression univariate FOOOF-LENS OS

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER bbp_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER slope_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER offset_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


* IDHmut, non-codeleted

USE ALL.
COMPUTE filter_$=(IDH_1p19q =1 AND Incl_tumor_probability =1).
VARIABLE LABELS filter_$ 'IDH_1p19q =0 AND Incl_tumor_probability =1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER bbp_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER slope_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER offset_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* Cox regression univariate FOOOF-LENS OS

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER bbp_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER slope_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER offset_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


* IDHwt

USE ALL.
COMPUTE filter_$=(IDH_1p19q =2 AND Incl_tumor_probability =1).
VARIABLE LABELS filter_$ 'IDH_1p19q =0 AND Incl_tumor_probability =1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER bbp_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER slope_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER offset_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* Cox regression univariate FOOOF-LENS OS

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER bbp_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER slope_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER offset_per_tumor 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


***** relationship with tumor volume while adjusting for tumor lat per subgoup

FILTER OFF.
USE ALL.
EXECUTE.

SPLIT FILE OFF.

USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

SORT CASES  BY IDH_1p19q.
SPLIT FILE LAYERED BY IDH_1p19q.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumor_volume_cm3
  /METHOD=ENTER bbp_per_tumor side.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumor_volume_cm3
  /METHOD=ENTER offset_per_tumor side.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumor_volume_cm3
  /METHOD=ENTER slope_per_tumor side.

** diffs in slope/offset/bbp  LENS per tumor population

FILTER OFF.
USE ALL.
EXECUTE.

SPLIT FILE OFF.

USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (slope_per_tumor offset_per_tumor bbp_per_tumor) GROUP (IDH_1p19q) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.


***** relationship with KPS 

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (slope_per_tumor offset_per_tumor bbp_per_tumor) GROUP (KPS_split) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

* get medians/ IQR

SORT CASES  BY KPS_split.
SPLIT FILE LAYERED BY KPS_split.

FREQUENCIES VARIABLES=bbp_per_tumor   offset_per_tumor  slope_per_tumor
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.

SPLIT FILE OFF.

* per subgroup (only slope significant)

SORT CASES  BY IDH_1p19q.
SPLIT FILE LAYERED BY IDH_1p19q.

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (slope_per_tumor) GROUP (KPS_split) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.


SORT CASES  BY KPS_split IDH_1p19q.
SPLIT FILE LAYERED BY KPS_split IDH_1p19q.

FREQUENCIES VARIABLES= slope_per_tumor
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.


* test whether AEDs / epi make a difference

SPLIT FILE OFF.

USE ALL.
COMPUTE filter_$=(Incl_tumor_probability=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (slope_per_tumor offset_per_tumor bbp_per_tumor) GROUP (epilepsy_dich) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

* get medians

SORT CASES  BY epilepsy_dich.
SPLIT FILE LAYERED BY epilepsy_dich.

FREQUENCIES VARIABLES= slope_per_tumor offset_per_tumor bbp_per_tumor
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.

SORT CASES  BY IDH_1p19q.
SPLIT FILE LAYERED BY IDH_1p19q.

FREQUENCIES VARIABLES= epilepsy_dich
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.
