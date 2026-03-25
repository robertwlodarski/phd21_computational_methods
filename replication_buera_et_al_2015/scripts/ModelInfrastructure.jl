# Content 
# 1. Parameters (structure and constructor)
# 2. Endogenous variables (structure and constructor)

# 1. Parameters (structure)
@with_kw struct ModelParameters
    # A. Parameters 
    σ::Float64          = 1.5           # Coefficient of RRA 
    δ::Float64          = 0.06          # Depreciation 
    γ::Float64          = 0.667         # Matching function parameter 
    α::Float64          = 0.33 * 0.79   # Capital share 
    θ::Float64          = 0.79 - α      # Labour share 
    η::Float64          = 5.25          # Pareto parameter 
    ψ::Float64          = 0.89          # Prob. of retaining productivity
    β::Float64          = 0.93          # Discount rate 
    λ::Float64          = 7.5           # Collateral constraint
    A::Float64          = 1.0           # Steady state productivity 

    # C. Grid sizes  
    Nᶻ::Int             = 40            # Productivity grids (number) 
    Nᵃ::Int             = 50            # Wealth grids (number)
    Nˡ::Int             = 40            # Employment grid

    # D. Productivity grid 
    # Pareto: f(z) = 1 - z^(-η), z ≥ 1
    z̲::Float64          = 1.0 + 1e-6    # Minimum productivity 
    z̅::Float64          = 10.0          # Maximum productivity 
    z⃗::Vector{Float64}  = zeros(Nᶻ)     # Grid 
    μ⃗::Vector{Float64}  = zeros(Nᶻ)     # PDF of the distribution 

    # E. Assets grid 
    a⃗::Vector{Float64}  = zeros(Nᵃ)     # Assets grid 
    a̲::Float64          = 0.0           # Minimum assets 
    a̅::Float64          = 700.0         # Maximum assets
    θᵃ::Float64         = 3.0           # Curvature of the assets grid 

    # F. Employment grid 
    l̲::Float64          = 0.0           # Minimum employment 
    l̅::Float64          = 2000.0        # Maximum employment 
    θˡ::Float64         = 3.0           # Curvature of employment grid 
    l⃗::Vector{Float64}  = zeros(Nˡ)     # Employment grid 

    # G. VFI-related parameters
    πˢᶜᵃˡᵉ::Float64     = 0.5           # Scale for the initial value function guess
    δʳᵉᶠ::Float64       = 0.01          # Tolerance for refining the grid
    𝒾̄ᵛᶠⁱ::Int           = 2500          # Maximum VFI iterations 
end 

# 1. Parameters (compiler)
function fnSetUpParameters(params::ModelParameters = ModelParameters())

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, z̲, z̅, η, a̲, a̅, θᵃ, l̲, l̅, θˡ = params

    # B. Productivity process 
    μ(z)    = 1 - z^(-η)
    z_of_p  = p -> (1 - p)^(-1/η)
    z⃗       = [LinRange(z_of_p(0.633), z_of_p(0.998), 38); z_of_p(0.999); z_of_p(0.9995)]
    μ⃗       = [μ(z⃗[1]); diff(μ.(z⃗))] ./ μ(z⃗[end])

    # C. Assets grid 
    a⃗       = a̲ .+ (a̅ .- a̲) .* (range(0,1,length=Nᵃ)).^θᵃ

    # D. Employment grid 
    l⃗       = l̲ .+ (l̅ .- l̲) .* (range(0, 1, length=Nˡ)).^θˡ

    # E. Save results 
    return reconstruct(params;
        z⃗   = z⃗,
        μ⃗   = μ⃗,
        a⃗   = a⃗, 
        l⃗   = l⃗
    )
end 
UsedParameters = fnSetUpParameters()

# 2. Endogenous variables preallocation (structure)
@with_kw mutable struct EndogenousVariables

    # A. Key values 
    𝐕::Matrix{Float64}      # Value function 
    𝔼𝐕::Matrix{Float64}     # Expected value function 
    𝐕ᵂ::Matrix{Float64}     # Value of working 
    𝐕ᴱ::Matrix{Float64}     # Value of entrepreneurship
    Π::Matrix{Float64}      # Firm profit  

    # B. Policy functions & indicators 
    𝐨::Matrix{Bool}         # Occupational choice (true = entrepreneurship)
    𝐤::Matrix{Float64}      # Capital policy function 
    𝐚::Matrix{Float64}      # Next period's assets 
    𝐜::Matrix{Float64}      # Consumption policy function 
    𝐥::Matrix{Float64}      # Labour policy function 
    𝕀ᶜ::Matrix{Bool}        # Constraint indicator 
    
    # C. Aggregate measures 
    G::Array{Float64, 4}    # Distribution 
    Kᵈ::Float64             # Capital demand 
    Kˢ::Float64             # Capital supplied 
    Lᵈ::Float64             # Labour demand 
    Lˢ::Float64             # Labour supply 
    Gᵉ::Float64             # Government expenditure

    # D. Prices 
    rₜ::Float64             # Interest rate 
    wₜ::Float64             # Wages 
    τₜ::Float64             # Taxes 
end

# 2. Endogenous variables preallocation (constructor)
function fnSetUpEndo(params::ModelParameters)
    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ = params 

    # B. Key values 
    𝐕       = zeros(Nᶻ, Nᵃ)
    𝔼𝐕      = zeros(Nᶻ, Nᵃ)
    𝐕ᵂ      = zeros(Nᶻ, Nᵃ)
    𝐕ᴱ      = zeros(Nᶻ, Nᵃ)
    Π       = zeros(Nᶻ, Nᵃ)

    # B. Policy functions & indicators 
    𝐨       = fill(true, Nᶻ, Nᵃ)
    𝐤       = zeros(Nᶻ, Nᵃ)
    𝐚       = zeros(Nᶻ, Nᵃ)
    𝐜       = zeros(Nᶻ, Nᵃ)
    𝐥       = zeros(Nᶻ, Nᵃ)
    𝕀ᶜ      = fill(true, Nᶻ, Nᵃ)
    
    # C. Aggregate measures 
    G       = zeros(Nᶻ, Nᵃ, Nˡ,2) # Last = unemployment
    Kᵈ      = 0.0 
    Kˢ      = 0.0
    Lᵈ      = 0.0
    Lˢ      = 0.0
    Gᵉ      = 0.0

    # D. Prices 
    rₜ      = 0.0
    wₜ      = 0.0
    τₜ      = 0.0

    # E. Returning 
    return EndogenousVariables(
      𝐕     = 𝐕,
      𝔼𝐕    = 𝔼𝐕,  
      𝐕ᵂ    = 𝐕ᵂ,    
      𝐕ᴱ    = 𝐕ᴱ, 
      Π     = Π,
      𝐨     = 𝐨, 
      𝐤     = 𝐤,
      𝐚     = 𝐚,
      𝐜     = 𝐜,
      𝐥     = 𝐥,
      𝕀ᶜ    = 𝕀ᶜ,
      G     = G,
      Kᵈ    = Kᵈ,
      Kˢ    = Kˢ,
      Lᵈ    = Lᵈ,
      Lˢ    = Lˢ, 
      Gᵉ    = Gᵉ,
      rₜ    = rₜ,
      wₜ    = wₜ,
      τₜ    = τₜ
    )
end 
Endo        = fnSetUpEndo(UsedParameters) 