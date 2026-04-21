# PHD21 Computational methods
# Replication: Elsby and Michaels (2013) Marginal Jobs, Heterogeneous Firms, and Unemployment Flows.
# January 2026

## 1. Packages
using DelimitedFiles, JLD2, DataInterpolations
using Preferences, Parameters, Accessors, StaticArrays, Adapt, QuantEcon, Printf
using Base.Cartesian, LinearAlgebra, SparseArrays, LoopVectorization, Interpolations
using Distributions, Random, StatsBase, FastGaussQuadrature, Optim, Roots, Dierckx
using BenchmarkTools, AllocCheck, MAT, Plots
Threads.nthreads()

## 2. Load macros, infrastructure, and functions
include("scripts/Macros.jl")
include("functions/InnerFunctions.jl")
include("scripts/ModelInfrastructure.jl")
include("functions/MainFunctions.jl")
include("functions/AggShocksFunctions.jl")

## 3. Search for the steady state at p=1.0
@time fSteadyState!(UsedParameters,Endo,1.0)

## 4. Using RTM to solve the problem 
@time fnSolveAggregateRTM!(UsedParameters, Endo, Simu, Lee; warm = true)
@save "results/rtm_warm_start.jld2" warm_q = Lee.q⃗ warm_N = Lee.N⃗ warm_Pi = Lee.Π̃ warm_n = Lee.n⃗

