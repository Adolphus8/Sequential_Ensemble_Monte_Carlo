# Sequential Ensemble Monte Carlo sampler
This repository presents a collection of tutorials (written in MATLAB) which seeks to demonstrate the implementation of the Sequential Ensemble Monte Carlo (SEMC) based on the literature by [Lye et. al (2022)](). Currently, 2 sets tutorials are presented here aimed at allowing users who are new such sampler find their footing around its concept and implementation. The details to these tutorials are as follows:

## Tutorials:

### 1) Numerical Example:
This example presents a SDoF Damped Oscillator whose spring stiffness and damping coefficient undergoes degradation following a random process. The objective is to perform an On-line Bayesian Inference and prediction over these parameters.

In this example, the proposed Sequential Ensemble Monte Carlo (SEMC) sampler is compared against the standard [Sequential Monte Carlo (SMC)](https://www.jstor.org/stable/4140600) sampler on the basis of its computation time, its estimates, and its robustness in moderating the acceptance rates within optimal bounds.

Note: Run "example_DampedOscillator_Part1.m" followed by "example_DampedOscillator_Part2.m"

### 2) Experimenal Example:
The experimental examples presents a SDoF single-storey structure subjected to Coulomb friction damping. Details to the experimental set-up an be found in the literature by [Marino et. al (2020)](https://doi.org/10.1007/s11071-019-05443-2).

In this example, the proposed SEMC sampler is used to infer the time-varying Coulomb friction as well as the time-invariant natural frequency and measurement errors based on actual experimental data. In addition, the sampler would also be used to compute the evidence of each Markov model used to model the degradation of the Coulomb friction and quantify the most appropriate model.

Note: Run "example_SDOF_System_Coulomb_Friction_Part1.m" followed by "example_SDOF_System_Coulomb_Friction_Part2.m"

## Reference(s):
* A. Lye, L. Marino, A. Cicirello, and E. Patelli (2022). Sequential Ensemble Monte Carlo Sampler for On-line Bayesian inference of Time-varying Parameter in Engineering Applications. *ASCE-ASME Journal of Risk and Uncertainty in Engineering Systems Part B: Mechanical Engineering*, 1–13, doi: [10.1115/1.4056934](https://asmedigitalcollection.asme.org/risk/article-abstract/doi/10.1115/1.4056934/1159641/Sequential-Ensemble-Monte-Carlo-Sampler-for-On?redirectedFrom=fulltext)
* A. Lye, A. Cicirello, and E. Patelli (2022). On-line Bayesian Inference for Structural Health Monitoring using Sequential Ensemble Monte Carlo Sampler. *In Proceedings of the 13th International Conference on Structural Safety and Reliability, 1*.

## Author:
* Name: Adolphus Lye
* Contact: adolphus.lye@liverpool.ac.uk
* Affiliation: Insitute for Risk and Uncertainty, University of Liverpool
