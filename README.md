# Dynamic Programming for Response-Adaptive Clinical Trials
Find the paper at http://snasrol.people.clemson.edu

1-The dose-response approximation for KGC and NIW approaches is a first order construction described in Section 3 of the paper. The solution method, a special kind of one-step look-ahead policy known as konwledge gradient, is described in Section 4. The objective function in KGC and NIW (parallel-palmetto) codes is to minimize the variance of the target dose, whereas in the NIW-Adaptive Randomization code, the objective is to maximize the probability of correctly identifying the target dose.

  KGC: Prior on mean responses is given by a multivariate normal distribution. Assumes mean dose responses are correlated and imposes a Guassian covariance structure over the covariance matrix. The observation variance is given and fixed.

  NIW (parallel-palmetto): Prior on mean responses is given by a normal inverse wishart distribution. No given correlation structure is assumed. Observation variance is unknown and may be different for each dose. The code is written in parallel to run on Clemson's Palmetto Cluster.
  
  NIW-Adaptive Randomization: Prior on mean responses is given by a normal inverse wishart distribution. No given correlation structure is assumed. Observation variance is unknown and may be different for each dose.

2-The dose-response approximation for NDLM-FFBS approach is a second order normal dynamic linear model described in details in Online Supplement Section 2. The solution method is similar to the knowledge gradient method, however, the posterior evaluation requires an implementation of the forward filtering, backward sampling algorithm.
  
  NDLM-FFBS: Prior distribution on NDLM parameters is assumed to be a bivariate normal. The variance is known here, but it may also be given an inverse Gamma prior. 
  
3-The dose-response approximation for the MCMC(Hill Curve) approach is given by a Hill equation which is a pre-spceified model to approximate sigmoid curves. This formulation and solution to this approach are described in details in Online Supplement Section 4.

  MCMC(Hill curve): Prior distribution on model paramteres are assumed to follow a uniform distribution. Posterior evaluation of the parameters are estimated by Markov Chaine Monte Carlo simulation.
