# PHD21 Computational methods
# Replication: Elsby and Michaels (2013) Marginal Jobs, Heterogeneous Firms, and Unemployment Flows.
# January 2026

# 1. Packages
using DelimitedFiles, JLD2, DataInterpolations
using Preferences, Parameters, Accessors, StaticArrays, Adapt, QuantEcon, Printf
using Base.Cartesian, LinearAlgebra, SparseArrays, LoopVectorization, Interpolations
using Distributions, Random, StatsBase, FastGaussQuadrature, Optim, Roots, Dierckx
using BenchmarkTools, AllocCheck, MAT, Plots
Threads.nthreads()

# 2. Load macros, infrastructure, and functions
include("scripts/Macros.jl")
include("functions/InnerFunctions.jl")
include("scripts/ModelInfrastructure.jl")
include("functions/MainFunctions.jl")
include("functions/AggShocksFunctions.jl")
include("functions/Plots.jl")
include("functions/Tables.jl")

# 3. Search for the steady state at p=1.0
@time fSteadyState!(UsedParameters,Endo,1.0)
fnSteadyStateTable(Endo,UsedParameters)

# 4. Using RTM to solve the problem 
# Add a max iteration cap so it doesn't run forever
# and saves on exit
try
    @time fnSolveAggregateRTM!(UsedParameters, Endo, Simu, Lee; warm = false)
catch e
    @warn "RTM interrupted: $e"
finally
    @save "results/rtm_warm_start_last_night.jld2" warm_q=Lee.q⃗ warm_N=Lee.N⃗ warm_Pi=Lee.Π̃ warm_n=Lee.n⃗
    println("💾 Final save completed")
end
# @save "results/rtm_warm_start.jld2" warm_q = Lee.q⃗ warm_N = Lee.N⃗ warm_Pi = Lee.Π̃ warm_n = Lee.n⃗

# 5. Start plotting 
Simu.N⃗ .= Lee.N⃗
Simu.q⃗ .= Lee.q⃗
Simu.f⃗ .= fUpdatedJobFindingRate.(Lee.q⃗, Ref(UsedParameters))
fnPlotBeveridgeCurve(Simu, UsedParameters)
fnPlotStateDepIRF(Simu, Lee, UsedParameters)
 fnCyclicalElasticityTable(Simu, Lee, UsedParameters)
 fnPlotSimulatedPaths(Simu, Lee, UsedParameters)
Endo.q
 @show std(lee.q⃗ⁱᵐᵖ), std(lee.q⃗)
@show std(lee.N⃗ⁱᵐᵖ), std(lee.N⃗)
