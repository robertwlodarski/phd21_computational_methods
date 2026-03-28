# PHD21 Computational method  
# Replication: Anatomy of a credit crunch: From capital to labor markets
# Authors: Buera, Jaef, and Shin (2015)
# March 2026
## 1. Packages & load functions 
using Parameters, Dierckx, Optim, NLsolve, QuantEcon, Plots, LinearAlgebra
include("scripts/ModelInfrastructure.jl")
include("scripts/Functions.jl")

@time fnSolveSteadyState!(UsedParameters, Endo)
