# Dynamic-DoseFinding
Dynamic response-adaptive design of dose-finding clinical trials (Knowledge gradient method, Normal dynamic linear method)

KGI: Knowledge gradient-Independent
Assumes mean dose responses are independent of each other and thus the prior covariance matrix is diagonal

KGC: Knowledge gradient-Correlated
Assumes mean dose responses are correlated and imposes a Guassian covariance structure over the covariance matrix.

NDLM-FFBS:
Assumes a second order normal dynamic linear model (general linear model) structure to approximate the dose-response curve. Uses
forward filtering backward samling (FFBS) to update posterior estimates of the NDLM model paramteres.
