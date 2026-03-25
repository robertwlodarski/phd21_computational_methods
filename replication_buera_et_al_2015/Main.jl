# PHD21 Computational method  
# Replication: Anatomy of a credit crunch: From capital to labor markets
# Authors: Buera, Jaef, and Shin (2015)
# March 2026

# PHD13 Public finance
# Replication: Life cycle model 
# March 2026

## 1. Packages & load functions 
using Parameters, QuantEcon, LinearAlgebra, Roots, Printf, Plots, Distributions, StatsBase, Random, Dierckx 
using Statistics, Plots
using DataFrames, Distributions, Random, GLM, Optim
include("scripts/ModelInfrastructure.jl")
include("scripts/Functions.jl")