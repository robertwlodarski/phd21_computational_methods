# Model infrastructure

# Content
#1. Params (struct & constructor)
#2. Endogenous variables (struct & constructor)
#3. Simulated variables  (struct & constructor)
#4. Krussel-Smith variabes (struct & constructor)

# 1.Parameters 
@with_kw struct ModelParameters
    # A. Grid and productivity parameters
    σ::Float64      = 0.25          # x shocks volatility 
    ξ::Float64                      # Pareto shape parameter (calculated in constructor)
    x̲::Float64                      # Pareto lower bound
    x̅::Float64                      # Pareto upper bound
    μₓ::Float64                     # Mean
    λ::Float64      = 0.043         # Shock arrival rate 
    p̄ₓ::Float64                     # Maximum probability
    
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
    
    # C. VFI-related parameters
    πˢᶜᵃˡᵉ::Float64 = 0.01          # Scale for the initial value function guess
    δʳᵉᶠ::Float64   = 0.01          # Tolerance for refining the grid
    n̅ˢ::Int         = 15            # The maximum number of spline interpolations
    N₁::Int         = 11            # First, sparse, segment of the VFI grid, number of elements
    N₂::Int         = 175           # Second, super dense, segment of the VFI grid, number of elements 
    N₃::Int         = 101           # Third, medium density, segment of the VFI grid, number of elements 
    N₄::Int         = 16            # Fourth, sparse, segment of the VFI grid, number of elements 
    N̅₁::Int         = 51            # First segment of the final grid, number of elements 
    N̅₂::Int         = 75            # Second segment of the final grid, number of elements 

    # D. Steady state computation settings
    q̅::Float64      = 0.95          # Upper bound for the job filling rate 
    q̲::Float64      = 0.05          # Lower bound for the job filling rate
    
    # E. Aggregate parameters 
    ρₚ::Float64     = 0.9925                # Productivity persistence 
    σ̃ₚ::Float64     = 0.0275                # Unconditional variance 
    σₚ::Float64     = σ̃ₚ * sqrt(1-ρₚ^2)     # Productivity update variance 
    Nₚ::Int         = 11                    # Number of Rouwenhorst grids
    P::Matrix{Float64}                      # Transition probability 
    p⃗::Vector{Float64}                      # Productivity grid  
    Ñₙ::Int         = 7                     # Number of aggregate employment grids
    Ñₜ::Int         = 11                    # Number of aggregate theta grids 

end

# 1. The constructor 
function fnSetUpParameters(; σ=0.25, Nₓ=45, ρₚ = 0.9925, σ̃ₚ = 0.0275, Nₚ = 11)
    # A. Calibrate Pareto
    σ̂               = 1.04 * σ
    ξ               = 1 + sqrt(1 + (1 / σ̂)^2)
    x̲               = (ξ - 1) / ξ
    x̅               = x̲ * ((1 - 0.9995)^(-1 / ξ))
    μₓ              = ξ * x̲ / (ξ - 1)
    
    # B. Generate grid (Gauss-Legendre)
    nodes, weights  = gausslegendre(Nₓ)
    x⃗               = nodes .* (x̅ - x̲) / 2.0 .+ (x̅ + x̲) / 2.0
    w⃗ₓ              = weights .* (x̅ - x̲) / 2.0
    
    # C. Calculate PDF and combined weights
    p̄ₓ              = 1 - (x̲ / x̅)^ξ
    𝑓x⃗              = (1 / p̄ₓ) .* (ξ * (x̲^ξ)) ./ (x⃗ .^ (ξ + 1))
    W⃗ₓ              = w⃗ₓ .* 𝑓x⃗          # Final weights for expectations, allowing E[V] = dot(W⃗ₓ, V)

    # D. Aggregate productivity items 
    σₚ              = σ̃ₚ * sqrt(1-ρₚ^2)
    ℳ𝒞              = rouwenhorst(Nₚ,ρₚ,σₚ) # Construct the Markov chain 
    P               = ℳ𝒞.p                  
    p⃗               = exp.(ℳ𝒞.state_values)      

    # E. Return the struct
    return ModelParameters(
        σ=σ, ξ=ξ, x̲=x̲, x̅=x̅, μₓ=μₓ, 
        Nₓ=Nₓ, x⃗=x⃗, W⃗ₓ=W⃗ₓ, p̄ₓ=p̄ₓ, P=P, p⃗=p⃗
    )
end

# 1. Run the parameters 
UsedParameters = fnSetUpParameters()

# 2. Endogenous variables preallocation
@with_kw mutable struct EndogenousVariables

    # A. Aggregate labour market values 
    θ::Float64      = 0         # Labour market tightness
    f::Float64      = 0         # Workers' contact rate 
    q::Float64      = 1         # Firms' contact rate 
    U::Float64      = 0         # Mass of unemployed 
    E::Float64      = 1         # Mass of employed 

    # B. Value functions 
    J::Matrix{Float64}          # Value function of a marginal job 
    Π::Matrix{Float64}          # Total firm value function
    Πᶜ::Matrix{Float64}         # Continuation value function
    Πᶠˡᵒʷ ::Matrix{Float64}     # Flow profit 
    Υ::Float64                  # Unemployment flow value 
    𝔼Π::Matrix{Float64}         # Expected value of firm value function
    n⃗::Vector{Float64}          # Employment grid
    R⃗::Vector{Float64}          # Firing threshold 
    ∂R⃗::Vector{Float64}         # Its partial derivative 
    R⃗ᵥ::Vector{Float64}         # Hiring threshold 
    ∂R⃗ᵥ::Vector{Float64}        # Its partial derivative

    # # C. Distributions 
    # 𝐆R⃗::Vector{Float64}         # Distribtion of firing threshold
    # 𝐆R⃗ᵥ::Vector{Float64}        # Distribtion of hiring threshold
    # 𝐇n⃗::Vector{Float64}         # Distribution of employment policy 
end

# 2. Constructor for endogenous variables 
function fnSetUpEndo(params::ModelParameters)
    # A. Unpack parameters 
    @unpack q̅ = params 

    # B. Initial guess for dimensions and labour grid  
    p⁰      = 1.0
    q⁰      = q̅
    f⁰      = fUpdatedJobFindingRate(q⁰,params)
    n⃗⁰      = fn⃗(params,p⁰,f⁰,q⁰)
    Nₓ      = params.Nₓ
    Nₙ      = length(n⃗⁰)

    # C. Allocate & return the structure 
    # Initialise with zeros, as the VFI will take care of the business 
    return EndogenousVariables(
        θ   = 0.0,
        f   = f⁰,
        q   = q⁰,
        U   = 0.0,
        E   = 0.0,
        Υ   = 0.0,
        n⃗   = n⃗⁰,
        J   = zeros(Nₓ,Nₙ),
        Π   = zeros(Nₓ,Nₙ),
        Πᶜ  = zeros(Nₓ,Nₙ),
        Πᶠˡᵒʷ   = zeros(Nₓ,Nₙ),
        𝔼Π  = zeros(Nₓ,Nₙ),
        R⃗   = zeros(Nₓ),
        ∂R⃗  = zeros(Nₓ),
        R⃗ᵥ  = zeros(Nₓ),
        ∂R⃗ᵥ = zeros(Nₓ)
    )
end 
Endo    = fnSetUpEndo(UsedParameters)

# 3. Simulated variables (structure)
# Important note:
# 𝒩⃗: caligraphic letters - state space
# N⃗: standard letters - simulation 
# Π̃: greek letters with a tilde - value functions with aggregate uncertainty

@with_kw mutable struct SimulationVariables

    # 1. Exogenous shocks
    p⃗̂::Vector{Float64}          # Simulated productivity vector 
    p⃗̂ᵢ::Vector{Int}             # Simulated productivity indices 
    
    # 2. Aggregate levels 
    N⃗::Vector{Float64}          # Vector of aggregate employment levels 
    Θ⃗::Vector{Float64}          # Vector of labour tightness values 
    q⃗::Vector{Float64}          # Vector of job-filling rates 
    f⃗::Vector{Float64}          # Vector of job-finding rates 

    # 3. Other flows 
    S⃗::Vector{Float64}          # Separations 
    M⃗::Vector{Float64}          # Matches 
    Y⃗::Vector{Float64}          # Output 
end 

# 3. Simulated variables (constructor)
function fnSetUpSimulated(UsedParameters; T = 52*75 +52*5, seed = 1997)

    # A. Unpacking and caring about reproducibility 
    @unpack P, p⃗, Nₚ            = UsedParameters
    Random.seed!(seed)

    # B. Productivity process 
    p⃗̂ᵢ              = zeros(Int,T)
    p⃗̂               = zeros(Float64,T)
    cs              = (Nₚ + 1) ÷ 2
    p⃗̂ᵢ[1]           = cs 
    p⃗̂[1]            = p⃗[cs]
    𝐹P              = cumsum(P,dims=2)
    for t in 1:(T-1)
        ns          = searchsortedfirst(𝐹P[cs,:],rand())
        p⃗̂ᵢ[t+1]     = ns 
        p⃗̂[t+1]      = p⃗[ns]
        cs          = ns 
    end 

    # C. Initialise empty vetors for aggregate endo variables 
    N⃗               = zeros(T) 
    Θ⃗               = zeros(T)
    q⃗               = zeros(T)
    f⃗               = zeros(T)
    S⃗               = zeros(T)
    M⃗               = zeros(T)
    Y⃗               = zeros(T)

    # Returning
    return SimulationVariables(
        p⃗̂ᵢ  = p⃗̂ᵢ,
        p⃗̂   = p⃗̂,
        N⃗   = N⃗,
        Θ⃗   = Θ⃗,
        q⃗   = q⃗,
        f⃗   = f⃗,
        S⃗   = S⃗,
        M⃗   = M⃗,
        Y⃗   = Y⃗
    )
end 

# 3. Simulate 
Simu    = fnSetUpSimulated(UsedParameters)   

#4. Krussel-Smith variabes (structure)
# Important note:
# 𝒩⃗: caligraphic latin letters - state space
# N⃗: standard latin letters - simulation  
# Π̃: greek letters with a tilde - value functions with aggregate uncertainty
@with_kw mutable struct KrussellSmithVariables

    # A. Regression parameters
    # Employment forecasting  
    ν₀::Float64     = 0.17634
    νₙ::Float64     = 0.92979
    νₚ::Float64     = 0.05397
    Rₙ::Float64     = 0.0
    # Tightness forecasting 
    θ₀::Float64     = -6.380265
    θₙ::Float64     = 1.338166
    θₚ::Float64     = 2.664425
    Rₜ::Float64     = 0.0
    # Collapsed forecasting parameters [updated using special functions]
    𝔼θₙ::Float64   
    𝔼θₚ::Float64    

    # B. Aggregate employment and tightness grids 
    𝒩⃗::Vector{Float64}
    ϴ⃗::Vector{Float64}
    𝓃⃗::Vector{Float64}

    # C. Multidimensional functions 
    # 5-dimensional: (x, p, N, Θ, n)
    # 4-dimensional: (x, p, N, n)
    Π̃::Array{Float64,5}             # Value function 
    𝔼Π̃::Array{Float64,4}            # Expected value function 
    Πᶜ::Array{Float64,5}            # Continuation value function 
end 

#4. Krussel-Smith variabes (constructor)
# Important note:
# 𝒩⃗: caligraphic latin letters - state space
# N⃗: standard latin letters - simulation  
# Π̃: greek letters with a tilde - value functions with aggregate uncertainty
function fnSetUpKS(UsedParameters)
end 