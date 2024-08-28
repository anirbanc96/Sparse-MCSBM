# Codes for experiments in [*Detecting Planted Partition in Sparse Multi-Layer Networks*](https://doi.org/10.1093/imaiai/iaae019)
This Github repository contains codes for experiments presented in our paper [Detecting Planted Partition in Sparse Multi-Layer Networks](https://doi.org/10.1093/imaiai/iaae019). We provide codes for generating figures $2-5$ from Section 5. The images, results, and corresponding generating scripts are organized as follows:

* Experiments from Section 5.1 (Dependence on Sample Size): [Over N](https://github.com/anirbanc96/Sparse-MCSBM/tree/main/Over%20N).
  * The results are obtained by running `Body.R` and images are generated by running `Plot.R`.
* Experiments from Section 5.2 (Dependence on Number of Networks): [Over M](https://github.com/anirbanc96/Sparse-MCSBM/tree/main/Over%20M).
  * The results are obtained by running `Body.R` and images are generated by running `Plot.R`.
* Experiments from Section 5.2 (Dependence on the ratio of the signal from covariates to the signal from the networks): [OverRatio](https://github.com/anirbanc96/Sparse-MCSBM/tree/main/OverRatio).
  * The results are obtained by running `Body.R` and images are generated by running `RatioPlot.R`.
* Experiments from Section 5.1 (Effect of data integration): [MethodsComparison](https://github.com/anirbanc96/Sparse-MCSBM/tree/main/MethodsComparison).
  * The results from Belief Propagation are obtained by running `BodyComp.R` from [MethodsComparison/BP](https://github.com/anirbanc96/Sparse-MCSBM/tree/main/MethodsComparison/BP)
  * The results from AMP are obtained by running `AMP.R` from [MethodsComparison/AMP](https://github.com/anirbanc96/Sparse-MCSBM/tree/main/MethodsComparison/AMP)
  * The results from DCMASE are obtained by running `Overlap.R` from [MethodsComparison/DCMASE](https://github.com/anirbanc96/Sparse-MCSBM/tree/main/MethodsComparison/DCMASE). We use the `R` implementation of DCMASE from [Joint Spectral Clustering in Multilayer Degree-Corrected Stochastic Blockmodels](https://arxiv.org/abs/2212.05053) given in the Github repository [jesusdaniel/dcmase/R](https://github.com/jesusdaniel/dcmase/tree/main/R).
