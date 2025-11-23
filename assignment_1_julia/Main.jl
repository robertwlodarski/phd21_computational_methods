# PHD21 Computational methods: Assignment 1
# Rob Włodarski 

## 1.       Load packages
using Plots
using Distributions
using Integrals
using Optim
using NLsolve
using Roots
using LinearAlgebra
using Infiltrator
using QuadGK

## 2.       Load functions 
include("_functions/InnerFunctions.jl")
include("_functions/Functions.jl")

## 3.       Set parameters
#           Exogeneous
const a     = 1.0
const α     = 0.3
const τ     = 0.15
const z̄     = 1.0
const A     = 1.0
const r     = 0.04
const β     = 0.96
#           Temporary
const η     = 6.4986
const χ     = 0.7744
const σ     = 0.3136
const b     = 0.0633

## 4.       Solve for equilibrium
@time wₑ,Tₑ       = fnWageSolver(A,α,r,z̄,a,b,τ,η,χ,β,σ)
