* Encoding: UTF-8.
* create pat chars for table 1

SORT CASES  BY group.
SPLIT FILE LAYERED BY group.

FREQUENCIES VARIABLES=sex kps_split lat IDH_1p19q epi_dich died progression
  /FORMAT=NOTABLE
  /ORDER=ANALYSIS.

DESCRIPTIVES VARIABLES=age
  /STATISTICS=MEAN STDDEV.

FREQUENCIES VARIABLES=tumorvol
  /FORMAT=NOTABLE
  /NTILES=4
  /ORDER=ANALYSIS.

USE ALL.
COMPUTE filter_$=(progression=1).
VARIABLE LABELS filter_$ 'progression=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES=PFS
  /FORMAT=NOTABLE
   /NTILES=4
  /ORDER=ANALYSIS.

USE ALL.
COMPUTE filter_$=(died=1).
VARIABLE LABELS filter_$ 'progression=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES=OS
  /FORMAT=NOTABLE
   /NTILES=4
  /ORDER=ANALYSIS.

FILTER OFF.
USE ALL.
EXECUTE.

*** output characteristics per mol group and cohort

USE ALL.
COMPUTE filter_$=(subgroup=1 OR subgroup = 2 OR subgroup = 3).
VARIABLE LABELS filter_$ 'progression=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

SORT CASES  BY group subgroup.
SPLIT FILE LAYERED BY group subgroup.

FREQUENCIES VARIABLES=sex kps_split lat IDH_1p19q epi_dich died progression
  /FORMAT=NOTABLE
  /ORDER=ANALYSIS.

DESCRIPTIVES VARIABLES=age
  /STATISTICS=MEAN STDDEV.

FREQUENCIES VARIABLES=tumorvol
  /FORMAT=NOTABLE
  /NTILES=4
  /ORDER=ANALYSIS.

USE ALL.
COMPUTE filter_$=(progression=1).
VARIABLE LABELS filter_$ 'progression=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES=PFS
  /FORMAT=NOTABLE
   /NTILES=4
  /ORDER=ANALYSIS.

USE ALL.
COMPUTE filter_$=(died=1).
VARIABLE LABELS filter_$ 'progression=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

FREQUENCIES VARIABLES=OS
  /FORMAT=NOTABLE
   /NTILES=4
  /ORDER=ANALYSIS.

FILTER OFF.
USE ALL.
EXECUTE.

SPLIT FILE OFF.

**** do molecular subgroup LENS analyses
    
USE ALL.
COMPUTE filter_$=(subgroup=1 OR subgroup = 2 OR subgroup = 3).
VARIABLE LABELS filter_$ 'progression=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.
    
NPTESTS 
  /INDEPENDENT TEST (LENSbbp LENSslope LENSoffset) GROUP (subgroup) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

* output medians/ranges

SORT CASES  BY subgroup.
SPLIT FILE LAYERED BY subgroup.

FREQUENCIES VARIABLES=LENSbbp LENSslope LENSoffset
  /FORMAT=NOTABLE
  /NTILES=4
  /ORDER=ANALYSIS.


*** do KPS LENS analyses

SPLIT FILE OFF.

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (LENSbbp LENSslope LENSoffset) GROUP (kps_split) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

SORT CASES  BY kps_split.
SPLIT FILE LAYERED BY kps_split.

FREQUENCIES VARIABLES=LENSbbp LENSslope LENSoffset
  /FORMAT=NOTABLE
  /NTILES=4
  /ORDER=ANALYSIS.

* significant for slope, so do posthoc per molecular subgroup

SORT CASES  BY subgroup.
SPLIT FILE LAYERED BY subgroup.

*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (LENSslope) GROUP (kps_split) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.


SORT CASES  BY subgroup kps_split.
SPLIT FILE LAYERED BY subgroup kps_split.

FREQUENCIES VARIABLES= LENSslope 
  /FORMAT=NOTABLE
  /NTILES=4
  /ORDER=ANALYSIS.

SPLIT FILE OFF.


*** Do tumor volume LENS analyses
    
SORT CASES  BY subgroup.
SPLIT FILE LAYERED BY subgroup.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumorvol
  /METHOD=ENTER LENSbbp lat.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumorvol
  /METHOD=ENTER LENSoffset lat.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT tumorvol
  /METHOD=ENTER LENSslope lat.

SPLIT FILE OFF.

*** Do epi LENS analyses
    
*Nonparametric Tests: Independent Samples. 
NPTESTS 
  /INDEPENDENT TEST (LENSbbp LENSslope LENSoffset) GROUP (epi_dich) 
  /MISSING SCOPE=ANALYSIS USERMISSING=EXCLUDE
  /CRITERIA ALPHA=0.05  CILEVEL=95.

SORT CASES  BY epi_dich.
SPLIT FILE LAYERED BY epi_dich.

FREQUENCIES VARIABLES= LENSbbp LENSslope LENSoffset 
  /FORMAT=NOTABLE
  /NTILES=4
  /ORDER=ANALYSIS.

SPLIT FILE OFF.

*** Do survival LENS analyses
    
SORT CASES  BY subgroup.
SPLIT FILE LAYERED BY subgroup.

COXREG PFS
  /STATUS=progression(1)
  /METHOD=ENTER ZLENSbbp 
   /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progression(1)
  /METHOD=ENTER ZLENSslope
   /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progression(1)
  /METHOD=ENTER ZLENSoffset
   /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER ZLENSbbp 
   /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER ZLENSslope
   /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER ZLENSoffset
   /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

