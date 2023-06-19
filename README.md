# RepMix
Repulsion, Chaos and Equilibrium in Mixture Models

# Authors
> Andrea Cremaschi, Singapore Institute for Clinical Sciences (SICS), Agency for Science, Technology and Research (A*STAR)

> Tim M. Wertz, Department of Mathematics, National University of Singapore (NUS)

> Maria De Iorio, Department of Paediatrics, Yong Loo Lin School of Medicine, National University of Singapore (NUS)

# Description
Mixture models are commonly used in applications with heterogeneity and overdispersion in the population, as they allow the identification of subpopulations. In the Bayesian framework, this entails the specification of suitable prior distributions for the weights and location parameters of the mixture. Widely used are Bayesian semi-parametric models based on mixtures with infinite or random number of components, such as Dirichlet process mixtures (Lo, 1984) or mixtures with random number of components (Miller and Harrison, 2018). Key in this context is the choice of the kernel for cluster identification. Despite their popularity, the flexibility of these models and prior distributions often does not translate into interpretability of the identified clusters. To overcome this issue, clustering methods based on repulsive mixtures have been recently proposed (Quinlan et al., 2021). The basic idea is to include a repulsive term in the prior distribution of the atoms of the mixture, which favours mixture locations far apart. This approach is increasingly popular and allows one to produce well-separated clusters, thus facilitating the interpretation of the results. However, the resulting models are usually not easy to handle due to the introduction of unknown normalising constants. Exploiting results from statistical mechanics, we propose in this work a novel class of repulsive prior distributions based on Gibbs measures. Specifically, we use Gibbs measures associated to joint distributions of eigenvalues of random matrices, which naturally possess a repulsive property. The proposed framework greatly simplifies the computations needed for the use of repulsive mixtures due to the availability of the normalising constant in closed form. 

This repository contains the code for implementing the proposed methodology. Specifically, we use the Air Quality dataset freely available in R.

# Contents
Simul_AirQuality.R
