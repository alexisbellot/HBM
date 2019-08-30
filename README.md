# Tree-based Bayesian Mixture model for Comepting Risks

This is an R implementation of the paper ["Tree-based Bayesian Mixture Model for Competing Risks"](http://proceedings.mlr.press/v84/bellot18a.html). 

In this project we  develop a semi-parametric Bayesian regression model for survival analysis with competing risks, which can be used for jointly assessing a patient’s risk of multiple (competing) adverse outcomes. The problem that motivated our approach is that in medical applications, therapies designed for patients at risk of multiple diseases need to account for the shared impact they may have on related diseases to ensure maximum overall well-being. Our algorithm is based on a Hierarchical Bayesian Mixture (HBM) model that describes survival paths in which a patient’s covariates influence both the estimation of the type of adverse event and the subsequent survival trajectory through Multivariate Random Forests. 

This repository includes the main model training function HBM.R, utility functions to simulate synthetic data and a self-contained demo application of our method. Currently only 2 competing risks are supported.

*Please cite the above paper if this resource is used in any publication*

## Requirements

* R version 3.5
* Packages: "randomForestSRC","flexsurv".

## First steps
To get started, check demo.R which will guide you through an application of our algorithm.

If you have questions or comments about anything regarding this work, please do not hesitate to contact [Alexis](https://alexisbellot.github.io/Website/)
