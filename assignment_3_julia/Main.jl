# PHD21 Computational methods: Assignment 1
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

## 2.       Load functions 
include("_functions/GridFunctions.jl")
include("_functions/InnerFunctions.jl")
include("_functions/Functions.jl")

## 3.       Set parameters and grids
#           Parameters
Base.@kwdef struct Parameters{T <: Real}
    α::T     = 0.36
    σ::T     = 0.20
    ρ::T     = 0.90
    β::T     = 0.96
    δ::T     = 0.08
    A::T     = 1.00
end 
p           = Parameters()
#           Grids
struct Grids
    vGridZ      ::Vector{Float64}
    mTransitionZ::Matrix{Float64}
    vGridA1     ::Vector{Float64}
    vGridA2     ::Vector{Float64}
end 
g           = Grids(fnTauchenLogNormal(p,3.0,7)..., fnGridMMV(0.0,150.0,50), fnGridMMV(0.0,150.0,100))
