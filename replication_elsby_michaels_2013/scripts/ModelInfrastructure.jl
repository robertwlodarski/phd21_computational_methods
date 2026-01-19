# Model infrastructure

# Content
#1. Params (struct & constructor)
#2. ExogVars (struct & constructor)


# 1. Parameters 
# Note: I follow Elsby and Michaels' notation from the paper (not replication package)
@with_kw mutable struct ModelParameters

    # A. Grid and productivity parameters 
    σ::Float64      = 0.25                          # x shocks volatility 
    σ̂::Float64      = 1.04*σ                        # Corrected volatility parameter
    ξ::Float64      = 1 + sqrt(1 + (1 / σ̂)^2)       # Pareto shape parameter
    x̲::Float64      = (ξ -1) / ξ                    # Pareto, lower bound
    x̅::Float64      = x̲ * ((1 - 0.9995)^(-1 / ξ))   # Pareto, upper bound
    μₓ::Float64     = ξ * x̲ / (ξ - 1)               # Mean of Pareto
    vₓ::Float64     = ξ * x̲ / ((ξ - 2) * (ξ - 1)^2) # Variance of Pareto 
    sₓ::Float64     = sqrt( vₓ )                    # St. dev. of Pareto 

    # B. Deep parameters
    ε::FLoat64      = 0.6                   # Matching elasticity
    μ::Float64      = 1.545 / 12            # Matching efficiency (CD form)
    r::Float64      = 0.012 / 12            # Discounting
    β::Float64      = 1 / (1 + r)           # Discounting
    α::Float64      = 0.60237               # Returns to scale
    η::Float64      = 0.443                 # Worker bargaining power for Stole & Zwiebel protocol
    L::Float64      = 3.5087                # Labour force size
    b::Float64      = 0.38735               # Flow value of unemployment 
    c::Float64      = 0.1327                # Vacancy posting cost
end