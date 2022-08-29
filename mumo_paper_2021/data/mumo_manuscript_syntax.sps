* Encoding: UTF-8.

USE ALL.

*** Main result-- multiple regression analysis single- and multilayer EC & EF

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA COLLIN TOL CHANGE
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN
  /DEPENDENT executive_functioning
  /METHOD=BACKWARD ec_meg_delta_fpn_mean ec_meg_theta_fpn_mean ec_meg_alpha1_fpn_mean
    ec_meg_alpha2_fpn_mean ec_meg_beta_fpn_mean ec_meg_gamma_fpn_mean ec_fmri_fpn_mean ec_dwi_fpn_mean
  /METHOD=ENTER multilayer_ec_fpn
  /SAVE ZRESID.

PPLOT
  /VARIABLES=ZRE_1
  /NOLOG
  /NOSTANDARDIZE
  /TYPE=Q-Q
  /FRACTION=BLOM
  /TIES=MEAN
  /DIST=NORMAL.



*** Main result-- multiple regression age and age squared & multilayer EC

COMPUTE age_sq=age * age.
EXECUTE.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA CHANGE
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN
  /DEPENDENT multilayer_ec_fpn
  /METHOD=ENTER age
  /METHOD=ENTER age_sq
  /SAVE ZRESID.

PPLOT
  /VARIABLES=ZRE_2
  /NOLOG
  /NOSTANDARDIZE
  /TYPE=Q-Q
  /FRACTION=BLOM
  /TIES=MEAN
  /DIST=NORMAL.

