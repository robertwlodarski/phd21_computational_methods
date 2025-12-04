# PhD21 Computational methods: Assignments

# Content

**Note**: 
1. The booklet is quite chaotic and very unfinished - I will fully update it *one* day. 
2. I've found a few errors in assignment 2, and I am updating the scripts to correct those. 

The repository contains a few folders. 
- `assignment_n`: In each of such folders, I store the `Matlab` code solving each of the assignments.
- `assignment_n_julia`: In each of such folders, I store the `Julia` version of the code. 
- `booklet`: This folder contains the written responses to each assignment.
- `replication_n`: In each of such folders, I store the code replicating one of the chosen papers.
- `extra_n`: In each of such folders, I store the code for an extra model that I replicated.
 
 Inside each `assignment_n` folder, I follow the same organisational structure:

```text
assignment_n                      	# Assignment folder
├── Main.m                  		# Master script
├── _scripts                      	# Code (used only once per assignment)
├── _functions                	    # Model implementations (e.g, equilibrium functions)
├── _plotting                       # Scripts used for preparing plots
├── _figures                  	    # Plots
```

# Assignments (non-graded)

1. RBC GE with labour participation.
2. A simple life-cycle PE model.
3. Aiyagari (QJE, 1994): Standard GE solution
4. Aiyagari (QJE, 1994): Transitional dynamics
5. RA model w/ aggregate uncertainty, Frisch labour supply, and adjustment costs
6. RA model w/ aggregate uncertainty, Frisch labour supply, and investment constraints
7. Krusell and Smith (1998) w/ aggregate uncertainty, Frisch labour supply, and adjustment costs

# Extras

On top of the required assignments, I've also used the course methods to solve the following models. 

1. Huggett (1993)
2. A standard search and matching DMP-style model.

# Replications (graded assignments)

1. Jarosch (2023): Unemployment trap model
2. Elsby and Michaels (2013): Marginal jobs and heterogeneous firms.

# Other items

## Requirements
No rocket science here. 
- `Matlab`: Apart from the standard tools, you need to download the *Global Optimisation Toolbox* for some of my calibration experiments.
- `LaTeX`: I use the standard VSCode for that. You may need to follow their instructions. 
- `Julia`: I try to translate the assignments to Julia, as well. I recommend using it with VSCode.

## Bibliography
- S. Rao Aiyagari, "Uninsured Idiosyncratic Risk and Aggregate Saving," *The Quarterly Journal of Economics*, Volume 109, Issue 3, August 1994, Pages 659–684, [https://doi.org/10.2307/2118417]
- Elsby, Michael W. L., and Ryan Michaels. 2013. "Marginal Jobs, Heterogeneous Firms, and Unemployment Flows." *American Economic Journal: Macroeconomics* 5 (1): 1–48.
- Jarosch, G. (2023), Searching for Job Security and the Consequences of Job Loss. *Econometrica*, 91: 903-942. [https://doi.org/10.3982/ECTA14008]
- Krusell, Per, and Anthony A. Smith, Jr. “Income and Wealth Heterogeneity in the Macroeconomy.” *Journal of Political Economy* 106, no. 5 (1998): 867–96. [https://doi.org/10.1086/250034]
- Lee, Hanbaek. 2025. "Global Nonlinear Solutions in Sequence Space and the Generalized Transition Function." *Working paper*
- Huggett, Mark, 1993. "The risk-free rate in heterogeneous-agent incomplete-insurance economies," *Journal of Economic Dynamics and Control*, Elsevier, vol. 17(5-6), pages 953-969.


