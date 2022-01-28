* Encoding: UTF-8.

*** get patient chars

FILTER OFF.
USE ALL.
EXECUTE.

SPLIT FILE OFF. 

FREQUENCIES VARIABLES=sex kps_split grade lat 
  /ORDER=ANALYSIS.


FREQUENCIES VARIABLES=age tumorvol
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.


USE ALL.
COMPUTE filter_$=(progressie = 1).
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
COMPUTE filter_$=(died = 1).
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

USE ALL.
COMPUTE filter_$=(inclusion_linda=1 OR inclusion_linda_idhmut=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.
    

SORT CASES  BY inclusion_linda.
SPLIT FILE LAYERED BY inclusion_linda.

FREQUENCIES VARIABLES=sex lat grade IDH_1p19q_interpolated KPS_split 
  /ORDER=ANALYSIS.

FREQUENCIES VARIABLES=AGE tumorvol 
  /NTILES=4
  /STATISTICS=STDDEV MEAN MEDIAN
  /ORDER=ANALYSIS.


USE ALL.
COMPUTE filter_$=((inclusion_linda=1 OR inclusion_linda_idhmut = 1) AND progressie = 1).
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
COMPUTE filter_$=((inclusion_linda=1 OR inclusion_linda_idhmut=1) AND died = 1).
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

**** diffs with HCs

T-TEST GROUPS=pt_vs_cont(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=age
  /CRITERIA=CI(.95).

CROSSTABS
  /TABLES=pt_vs_cont BY sex
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT
  /COUNT ROUND CELL
  /METHOD=EXACT TIMER(5).


* Cox regression univariate LENS-bbp & FOOOF PFS

* IDH mut

USE ALL.
COMPUTE filter_$=(inclusion_linda_idhmut=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
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

* None of them died, so no OS
* no KM because non-significant

***** relationship with tumor volume

FILTER OFF.
USE ALL.
EXECUTE.

SORT CASES  BY inclusion_linda.
SPLIT FILE LAYERED BY inclusion_linda.

* regress while adjusting for tumor lateralization

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumorvol
  /METHOD=ENTER lat bbp_per_tumor.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumorvol
  /METHOD=ENTER lat offset_per_tumor.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumorvol
  /METHOD=ENTER lat slope_per_tumor.



* difference between groups? 

FILTER OFF.
USE ALL.
EXECUTE.

SPLIT FILE OFF.

USE ALL.
COMPUTE filter_$=(inclusion_linda = 1 OR inclusion_linda = 2).
VARIABLE LABELS filter_$ 'IDH_1p19q_interpolated = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (slope_per_tumor offset_per_tumor bbp_per_tumor) GROUP (inclusion_linda) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.



***** relationship with KPS 

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (slope_per_tumor offset_per_tumor bbp_per_tumor) GROUP (KPS_split) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

SORT CASES  BY kps_split.
SPLIT FILE LAYERED BY kps_split.


FREQUENCIES VARIABLES=bbp_per_tumor offset_per_tumor slope_per_tumor 
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.


** epilepsy? 

NPTESTS 
  /INDEPENDENT TEST (slope_per_tumor offset_per_tumor bbp_per_tumor) GROUP (aed_dich) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

* get group-level medians

SORT CASES  BY aed_dich.
SPLIT FILE LAYERED BY aed_dich.

FREQUENCIES VARIABLES=bbp_per_tumor offset_per_tumor slope_per_tumor 
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.

* bbp is significant, so do subgroup analysis (posthoc)

SPLIT FILE OFF.

USE ALL.
COMPUTE filter_$=(inclusion_linda=1).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.
    
NPTESTS 
  /INDEPENDENT TEST (bbp_per_tumor) GROUP (aed_dich) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

SORT CASES  BY aed_dich.
SPLIT FILE LAYERED BY aed_dich.

FREQUENCIES VARIABLES=bbp_per_tumor 
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.

SPLIT FILE OFF. 

USE ALL.
COMPUTE filter_$=(inclusion_linda=2).
VARIABLE LABELS filter_$ 'Incl_tumor_probability=1 (FILTER)'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

NPTESTS 
  /INDEPENDENT TEST (bbp_per_tumor) GROUP (aed_dich) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

SORT CASES  BY aed_dich.
SPLIT FILE LAYERED BY aed_dich.

FREQUENCIES VARIABLES=bbp_per_tumor 
  /NTILES=4
  /STATISTICS= MEDIAN
  /ORDER=ANALYSIS.


