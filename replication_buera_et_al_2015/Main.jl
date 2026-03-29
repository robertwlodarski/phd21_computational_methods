# PHD21 Computational method  
# Replication: Anatomy of a credit crunch: From capital to labor markets
# Authors: Buera, Jaef, and Shin (2015)
# March 2026
# Remember to change the number of threads 
## 1. Packages & load functions 

using Parameters, Dierckx, Optim, NLsolve, QuantEcon, Plots, LinearAlgebra, Roots
Threads.nthreads()
include("scripts/ModelInfrastructure.jl")
include("scripts/Functions.jl")


@time fnSolveSteadyState!(UsedParameters, Endo)
