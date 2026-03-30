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
    Nᵃ::Int             = 70            # Wealth grids (number)
    Nˡ::Int             = 50            # Employment grid
    Nᵘ::Int             = 2             # Unemployment and other states grid 

    # D. Productivity grid 
    # Pareto: f(z) = 1 - z^(-η), z ≥ 1
    z̲::Float64          = 1.0 + 1e-6    # Minimum productivity 
    z̅::Float64          = 4.31          # Maximum productivity 
    z⃗::Vector{Float64}  = zeros(Nᶻ)     # Grid 
    μ⃗::Vector{Float64}  = zeros(Nᶻ)     # PDF of the distribution 

    # E. Assets grid 
    a⃗::Vector{Float64}  = zeros(Nᵃ)     # Assets grid 
    a̲::Float64          = 0.0           # Minimum assets 
    a̅::Float64          = 7500.0        # Maximum assets
    θᵃ::Float64         = 3.0           # Curvature of the assets grid
    c̲::Float64          = 1e-6          # "Zero" consumption  

    # F. Employment grid 
    l̲::Float64          = 0.0           # Minimum employment 
    l̅::Float64          = 5000.0        # Maximum employment 
    θˡ::Float64         = 3.0           # Curvature of employment grid 
    l⃗::Vector{Float64}  = zeros(Nˡ)     # Employment grid 

    # G. VFI-related and distribution-related parameters
    δᵛᶠⁱ::Float64       = 1e-2          # VFI iteration tolerance
    𝒾̄ᵛᶠⁱ::Int           = 2000          # Maximum VFI iterations 
    λᵛᶠⁱ::Float64       = 0.0           # Updating loading (VFI)
    δᵈⁱˢᵗ::Float64      = 1e-8          # Distribution iteration tolerance
    λᵈⁱˢᵗ::Float64      = 0.5           # Updating loading for distributions

    # H. GE bounds for the prices 
    w̲::Float64          = 0.99              # Minimum wage 
    w̅::Float64          = 2.50              # Maximum wage 
    r̲::Float64          = 0.005             # Minimum interest rate 
    r̅::Float64          = 1.0 / β - 1.02    # Maximum interest rate 
    τ̲::Float64          = 0.99 * 0.05       # Minimum tax 
    τ̅::Float64          = 2.50*0.25         # Maximum tax 

    # I. Labour market loop updating 
    δᴸ::Float64         = 5*1e-3            # Tolerance 
    κᴸ::Float64         = 0.40              # Updates
    λᵗ::Float64         = 0.50              # Tax update
    κᵗ::Float64         = 0.70              # Updates
    δᵗ::Float64         = 1e-3              # Tolerance
    δʳ::Float64         = 1e-4              # Tolerance
    κʳ::Float64         = 0.70              # Updates   
end 

# 1. Parameters (compiler)
function fnSetUpParameters(params::ModelParameters = ModelParameters())

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ,Nˡ, z̲, z̅, η,α, a̲, a̅, θᵃ, l̲, l̅, θˡ,θ = params

    # B. Productivity process 
    μ(z)            = 1 - z^(-η)
    z_of_p(p)       = (1 - p)^(-1/η)
    step_p          = 0.9995 / (Nᶻ - 1)
    p_bounds        = collect(range(step_p, 0.9995, length=Nᶻ-1))    
    z⃗               = zeros(Nᶻ)
    μ⃗               = zeros(Nᶻ)
    μ⃗[1]            = p_bounds[1]
    z⃗[1]            = z_of_p(p_bounds[1] / 2.0)
    for i in 2:(Nᶻ - 1)
        μ⃗[i]        = p_bounds[i] - p_bounds[i-1]
        p_mid       = (p_bounds[i-1] + p_bounds[i]) / 2.0
        z⃗[i]        = z_of_p(p_mid)
    end
    μ⃗[Nᶻ]           = 1.0 - p_bounds[end]
    ν               = α + θ                 
    ξ               = 1.0 / (1.0 - ν)       
    if η <= ξ
        error("Pareto tail integral diverges.")
    end
    tail_multiplier = (η / (η - ξ)) ^ (1.0 / ξ)
    z⃗[Nᶻ]           = z_of_p(p_bounds[end]) * tail_multiplier

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
    𝐚ʷ::Matrix{Float64}      # → of worker
    𝐚ᵉ::Matrix{Float64}      # → of entrepreneur 
    𝐜::Matrix{Float64}      # Consumption policy function 
    𝐜ʷ::Matrix{Float64}      # → of worker 
    𝐜ᵉ::Matrix{Float64}      # → of entrepreneur 
    𝐥::Matrix{Float64}      # Labour policy function 
    𝕀ᶜ::Matrix{Bool}        # Constraint indicator 
    
    # C. Aggregate measures 
    g::Array{Float64, 4}    # Distirbution (PDF)
    g̃::Array{Float64, 3}    # Marginal distribution (PDF)
    JD::Float64             # Job destruction 
    S::Float64              # Switchers from entrepreneurship to working 
    D::Float64              # Jobs reduced 
    Kᵈ::Float64             # Capital demand 
    Kˢ::Float64             # Capital supplied 
    Lᵈ::Float64             # Labour demand 
    Lˢ::Float64             # Labour supply 
    Gᵉ::Float64             # Government expenditure
    U::Float64              # Unemployed workers 
    M::Float64              # Matches 
    W::Float64              # Employed workers 
    E::Float64              # Entrepreneurs
    
    # D. Prices 
    rₜ::Float64             # Interest rate 
    wₜ::Float64             # Wages 
    τₜ::Float64             # Taxes 
end

# 2. Endogenous variables preallocation (constructor)
function fnSetUpEndo(params::ModelParameters)
    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ,Nᵘ = params 

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
    𝐚ᵉ      = zeros(Nᶻ, Nᵃ)
    𝐚ʷ      = zeros(Nᶻ, Nᵃ)
    𝐜       = zeros(Nᶻ, Nᵃ)
    𝐜ʷ      = zeros(Nᶻ, Nᵃ)
    𝐜ᵉ      = zeros(Nᶻ, Nᵃ)
    𝐥       = zeros(Nᶻ, Nᵃ)
    𝕀ᶜ      = fill(true, Nᶻ, Nᵃ)
    
    # C. Aggregate measures 
    g       = ones(Nᶻ, Nᵃ, Nˡ, Nᵘ) ./ (Nᶻ * Nᵃ * Nˡ * Nᵘ)
    g̃       = zeros(Nᶻ, Nᵃ, Nˡ)
    JD      = 0.0
    D       = 0.0
    S       = 0.0
    Kᵈ      = 0.0 
    Kˢ      = 0.0
    Lᵈ      = 0.0
    Lˢ      = 0.0
    Gᵉ      = 0.0
    U       = 0.0
    M       = 0.0
    W       = 0.0
    E       = 0.0

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
      𝐚ᵉ    = 𝐚ᵉ,
      𝐚ʷ    = 𝐚ʷ,
      𝐜     = 𝐜,
      𝐜ʷ    = 𝐜ʷ,
      𝐜ᵉ    = 𝐜ᵉ,
      𝐥     = 𝐥,
      𝕀ᶜ    = 𝕀ᶜ,
      g     = g,
      g̃     = g̃,
      JD    = JD,
      D     = D, 
      S     = S,
      Kᵈ    = Kᵈ,
      Kˢ    = Kˢ,
      Lᵈ    = Lᵈ,
      Lˢ    = Lˢ, 
      Gᵉ    = Gᵉ,
      U     = U, 
      M     = M,
      W     = W,
      E     = E,
      rₜ    = rₜ,
      wₜ    = wₜ,
      τₜ    = τₜ
    )
end 
Endo        = fnSetUpEndo(UsedParameters) 