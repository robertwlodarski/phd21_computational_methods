# PHD21 Computational methods
# Extra: DMP model 
# Rob Włodarski 

## 1.       Load packages
using Distributions
using Integrals
using Optim
using NLsolve
using Roots
using LinearAlgebra
using Infiltrator
using QuadGK
using Parameters
using Pkg
using SparseArrays
using Printf

## 2.       Load functions 
include("_functions/InnerFunctions.jl")
include("_functions/Functions.jl")

## 3.       Set parameters and grids
#           Parameters
Base.@kwdef struct ModelParameters{T <: Real}
    μ::T        = 0.2910        # Match efficiency
    ξ::T        = 0.6353        # Matching function parameters
    β::T        = 0.99^(1/3)    # Discount factor
    σ::T        = 1.0           # Intertemporal elasticity of substitution
    λ::T        = 0.0283        # Exogenous separations
    κ::T        = 0.0667        # Vacancy posting cost
    b::T        = 0.6720        # Unemployment benefit
    η::T        = 0.6353        # Bargaining share of workers
end 
p           = ModelParameters()

## 4.       Solve the steady state model
@time θ1, w1    = fnPlainDMPSteadyStateNaive(1.0,p)
@time θ2, w2    = fnPlainDMPSteadyStateRoot(1.0,p)