# Important note on notation with aggregate uncertainty 
# 𝒩⃗: caligraphic letters:           State space
# N⃗: standard letters:              Simulation  
# Π̃: greek letters with a tilde:    Value functions

# Content: 
# 1A. Collapsed forecasting parameters 
# 1B. Prepare aggregate N grid 

# 1A. Collapsed forecasting parameters 
function fnCollapsedForecastingParameters(ν₀,νₙ,νₚ,θ₀,θₙ,θₚ)

    # A. Compute the simplification 
    𝔼θₙ         = θₙ*νₙ
    𝔼θₚ         = θₚ+θₙ*νₚ
    𝔼θ₀         = θ₀+θₙ*ν₀

    # B. Return values
    return 𝔼θₙ, 𝔼θₚ, 𝔼θ₀
end 

# 1B. Prepare aggregate grids
function fnAggregateGridsStateSpace!(params,KS,ν₀,νₙ,νₚ,θ₀,θₙ,θₚ)

    # A. Unpacking business 
    @unpack p⃗,Ñₙ,Ñₜ = params 

    # B. Prepare the "collapsed" forecasting parameters
    KS.𝔼θₙ, KS.𝔼θₚ, KS.𝔼θ₀   = fnCollapsedForecastingParameters(ν₀,νₙ,νₚ,θ₀,θₙ,θₚ)

    # C. Aggregate employment state space 
    𝒩⃗̄               = max(3.333,(ν₀+νₚ*maximum(p⃗))/(1-νₙ))
    𝒩̲⃗               = min(3.220,(ν₀+νₚ*minimum(p⃗))/(1-νₙ))
    Δ𝒩              = (𝒩⃗̄-𝒩̲⃗)/(Ñₙ-1)
    KS.𝒩⃗            = collect(𝒩̲⃗:Δ𝒩:𝒩⃗̄)
end 
