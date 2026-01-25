# Model infrastructure

# Content
#1. Params (struct & constructor)
#3. ExogVars (struct & constructor)


using Parameters, FastGaussQuadrature, LinearAlgebra

# 1.Parameters 

@with_kw struct ModelParameters
    # A. Grid and productivity parameters
    σ::Float64      = 0.25          # x shocks volatility 
    ξ::Float64                      # Pareto shape parameter (calculated in constructor)
    x̲::Float64                      # Pareto lower bound
    x̅::Float64                      # Pareto upper bound
    μₓ::Float64                     # Mean
    λ::Float64      = 0.043         # Shock arrival rate 
    
    # Grid vectors
    Nₓ::Int         = 45            # Number of grids 
    x⃗::Vector{Float64}              # The productivity Grid
    W⃗ₓ::Vector{Float64}             # The combined integration weights (w * pdf)
    Δₙ::Float64     = 0.02          # Common step for the endogenously-set
    s̅ₙ::Float64     = 1.25          # Grid safety upper multiplier 
    s̲ₙ::Float64     = 5             # Grid safety bottom multiplier   

    # B. Deep parameters
    ε::Float64      = 0.6           # Matching elasticity
    μ::Float64      = 1.545 / 12    # Matching efficiency
    r::Float64      = 0.012 / 12    # Discounting
    β::Float64      = 1 / (1 + r)   # Discount factor
    α::Float64      = 0.60237       # Returns to scale
    η::Float64      = 0.443         # Bargaining power
    L::Float64      = 3.5087        # Labour force
    b::Float64      = 0.38735       # Unemployment flow
    c::Float64      = 0.1327        # Vacancy cost
    pₛₛ::Float64    = 1             # Steady state productivity 
end

# 1. The constructor 
function setup_parameters(; σ=0.25, Nₓ=45)
    # 1. Calibrate Pareto
    σ̂               = 1.04 * σ
    ξ               = 1 + sqrt(1 + (1 / σ̂)^2)
    x̲               = (ξ - 1) / ξ
    x̅               = x̲ * ((1 - 0.9995)^(-1 / ξ))
    μₓ              = ξ * x̲ / (ξ - 1)
    
    # 2. Generate grid (Gauss-Legendre)
    nodes, weights  = gausslegendre(Nₓ)
    x⃗               = nodes .* (x̅ - x̲) / 2.0 .+ (x̅ + x̲) / 2.0
    w⃗ₓ              = weights .* (x̅ - x̲) / 2.0
    
    # 3. Calculate PDF and combined weights
    p̄ₓ              = 1 - (x̲ / x̅)^ξ
    pdf_x           = (1 / p̄ₓ) .* (ξ * (x̲^ξ)) ./ (x⃗ .^ (ξ + 1))
    W⃗ₓ              = w⃗ₓ .* pdf_x        # Final weights for expectations, allowing E[V] = dot(W⃗ₓ, V)

    # 4. Return the struct
    return ModelParameters(
        σ=σ, ξ=ξ, x̲=x̲, x̅=x̅, μₓ=μₓ, 
        Nₓ=Nₓ, x⃗=x⃗, W⃗ₓ=W⃗ₓ
    )
end

# 1. Run the parameters 
UsedParameters = setup_parameters()

# 2. Endogenous variables preallocation

@with_kw mutable struct EndogenousVariables

    # A. Aggregate labour market values 
    θ::Float64      = 0         # Labour market tightness
    f::Float64      = 0         # Workers' contact rate 
    q::Float64      = 1         # Firms' contact rate 
    U::Float64      = 0         # Mass of unemployed 
    E::Float64      = 1         # Mass of employed 

    # B. Value functions 
    J::AbstractArray{Float64,2}     # Value function of a marginal job 
    Π::AbstractArray{Float64,2}     # Total firm value function 
    Υ::Float64                      # Unemployment flow value 

end