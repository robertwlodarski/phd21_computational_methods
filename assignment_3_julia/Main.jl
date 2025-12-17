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
using Parameters
using Pkg

## 2.       Load functions 
include("_functions/InnerFunctions.jl")
include("_functions/Functions.jl")

## 3.       Set parameters
Base.@kwdef struct Parameters{T <: Real}
    α::T     = 0.36
    σ::T     = 0.20
    ρ::T     = 0.90
    β::T     = 0.96
    δ::T     = 0.08
    A::T     = 1.00
end 
p           = Parameters()

