# PHD21 Computational method  
# Replication: Anatomy of a credit crunch: From capital to labor markets
# Authors: Buera, Jaef, and Shin (2015)
# March 2026
# Remember to change the number of threads 
## 1. Packages & load functions 
using Parameters, Dierckx, Optim, NLsolve, QuantEcon, Plots, LinearAlgebra, Roots, Printf, Interpolations
Threads.nthreads()
include("scripts/ModelInfrastructure.jl")
include("scripts/Functions.jl")
include("scripts/FunctionsMITShock.jl")
include("scripts/PlottingMIT.jl")
include("scripts/Tables.jl")

## 2. Solve the steady state model 
@time fnSolveSteadyState!(UsedParameters, Endo)
fnPrintCalibrationElements(UsedParameters, Endo)
fnPrintCalibrationLaTeX(UsedParameters, Endo; variant = :grid, path = joinpath("tables", "calibration.tex"))

## 2B. Run the model for the values of solved prices. I does bad. 
# Endo.wₜ = 1.532166
# Endo.τₜ = 0.101028
# Endo.rₜ = 0.009158
# fnVFI!(UsedParameters, Endo)
# fnAggregateStates!(UsedParameters, Endo)
# fnPrintCalibrationElements(UsedParameters, Endo)

## 3. Solve the MIT transition 
@time fnSolveMIT!(UsedParameters, EndoMIT, Endo, ConstantTechnology, CollateralShock)
plt = fnPlotMITResults(UsedParameters, EndoMIT, Endo,ConstantTechnology)
savefig(plt, "plots/MIT_transitions.pdf")
