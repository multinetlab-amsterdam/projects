* Encoding: UTF-8.

* Cox regression univariate LENS-bbp & FOOOF PFS

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER bbp_per_tumor IDH_1p19q tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER slope_per_tumor IDH_1p19q tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER offset_per_tumor IDH_1p19q tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* Cox regression univariate FOOOF-LENS OS

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER bbp_per_tumor IDH_1p19q tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER slope_per_tumor IDH_1p19q tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER offset_per_tumor IDH_1p19q tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


****** too many covars? do per subgroup
*idhmut

USE ALL.
COMPUTE filter_$=(IDH_1p19q_interpolated = 1).
VARIABLE LABELS filter_$ 'IDH_1p19q_interpolated = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER bbp_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER slope_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER offset_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* Cox regression univariate FOOOF-LENS OS

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER bbp_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER slope_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER offset_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* idhwt

USE ALL.
COMPUTE filter_$=(IDH_1p19q_interpolated = 2).
VARIABLE LABELS filter_$ 'IDH_1p19q_interpolated = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER bbp_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER slope_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER offset_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* Cox regression univariate FOOOF-LENS OS

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER bbp_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER slope_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG OS
  /STATUS=died(1)
  /METHOD=ENTER offset_per_tumor tumorvol age kps_split
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* median split

FILTER OFF.
USE ALL.
EXECUTE.

FREQUENCIES VARIABLES=slope_per_tumor offset_per_tumor bbp_per_tumor
  /NTILES=4
  /ORDER=ANALYSIS.

*LD: adjusted these numbers after combining AMS & BOS
*median_FOOOF_slope =1.1023
*median_FOOOF_offset=7.1466
* median_bbp = 31991.97

DO IF  (slope_per_tumor > 1.1023). 
COMPUTE slope_median_split=1.
ELSE IF  (slope_per_tumor <=1.1023).
COMPUTE  slope_median_split=0.
END IF.
EXECUTE.

DO IF  (offset_per_tumor > 7.1466). 
COMPUTE offset_median_split=1.
ELSE IF  (offset_per_tumor <=7.1466).
COMPUTE offset_median_split=0.
END IF.
EXECUTE.

DO IF  (bbp_per_tumor >31991.97). 
COMPUTE bbp_median_split=1.
ELSE IF  (bbp_per_tumor <=31991.97).
COMPUTE bbp_median_split=0.
END IF.
EXECUTE.

* Define Variable Properties.
*slope_median_split.
VALUE LABELS slope_median_split
  .00 'low slope tumor location'
  1.00 'high slope tumor location'.
*offset_median_split.
VALUE LABELS offset_median_split
  .00 'low offset tumor location'
  1.00 'high offset tumor location'.
*bbp_median_split.
VALUE LABELS bbp_median_split
  .00 'low broadband power tumor location'
  1.00 'high broadband power tumor location'.
EXECUTE.


* Create Kaplan Meier  plots

KM PFS BY bbp_median_split
  /STATUS=progressie(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL.

KM PFS BY slope_median_split
  /STATUS=progressie(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL.

KM PFS BY offset_median_split
  /STATUS=progressie(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL.

KM OS BY bbp_median_split
  /STATUS=died(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL.

KM OS BY slope_median_split
  /STATUS=died(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL.

KM OS BY offset_median_split
  /STATUS=died(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL.


* Check PFS with bbp with covars per subgroup

* idhmut non-codel

USE ALL.
COMPUTE filter_$=(IDH_1p19q_interpolated = 1).
VARIABLE LABELS filter_$ 'IDH_1p19q_interpolated = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER age tumorvol bbp_per_tumor
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER age tumorvol slope_per_tumor
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER age tumorvol offset_per_tumor
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

* idhwt

USE ALL.
COMPUTE filter_$=(IDH_1p19q_interpolated = 2).
VARIABLE LABELS filter_$ 'IDH_1p19q_interpolated = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER age tumorvol kps_split bbp_per_tumor
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER age tumorvol kps_split slope_per_tumor
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG PFS
  /STATUS=progressie(1)
  /METHOD=ENTER age tumorvol kps_split offset_per_tumor
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

KM PFS BY bbp_median_split
  /STATUS=progressie(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL.



***** relationship with tumor volume

FILTER OFF.
USE ALL.
EXECUTE.

SORT CASES  BY IDH_1p19q.
SPLIT FILE LAYERED BY IDH_1p19q.

NONPAR CORR
  /VARIABLES=tumorvol slope_per_tumor offset_per_tumor bbp_per_tumor
  /PRINT=SPEARMAN TWOTAIL NOSIG
  /MISSING=PAIRWISE.

* difference between groups? 

FILTER OFF.
USE ALL.
EXECUTE.

SPLIT FILE OFF.


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

