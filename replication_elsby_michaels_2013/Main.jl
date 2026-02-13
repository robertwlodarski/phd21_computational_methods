# PHD21 Computational methods
# Replication: Elsby and Michaels (2013) Marginal Jobs, Heterogeneous Firms, and Unemployment Flows.
# January 2026

## 1. Packages
using Preferences, Parameters, Accessors, StaticArrays, Adapt
using Base.Cartesian, LinearAlgebra, SparseArrays, LoopVectorization, Interpolations
using Distributions, Random, StatsBase, FastGaussQuadrature, Optim, Roots, Dierckx
using BenchmarkTools, AllocCheck, MAT

## 2. Load macros, infrastructure, and functions
include("scripts/Macros.jl")
include("functions/InnerFunctions.jl")
include("scripts/ModelInfrastructure.jl")
include("functions/MainFunctions.jl")

## 3. Search for the steady state at p=1.0
q̂,f̂,N̂,Ŷ,Ŝ,M̂,Â =@time fSteadyState(UsedParameters,Endo,1.0)
