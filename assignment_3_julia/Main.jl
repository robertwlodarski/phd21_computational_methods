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
const α     = 0.36
const σ     = 0.20
const ρ     = 0.90
const β     = 0.96
const δ     = 0.08
const A     = 1.00


