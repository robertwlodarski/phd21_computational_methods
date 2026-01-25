# PHD21 Computational methods
# Replication: Elsby and Michaels (2013) Marginal Jobs, Heterogeneous Firms, and Unemployment Flows.
# January 2026

## 1. Packages
Pkg.status()
using Preferences, Parameters, Accessors, StaticArrays, Adapt
using Base.Cartesian, LinearAlgebra, SparseArrays, LoopVectorization
using Distributions, Random, StatsBase, FastGaussQuadrature
using BenchmarkTools, AllocCheck, MAT 

## 2. Load macros, infrastructure, and functions
include("scripts/Macros.jl")
include("scripts/ModelInfrastructure.jl")
include("_functions/InnerFunctions.jl")
include("_functions/Functions.jl")


