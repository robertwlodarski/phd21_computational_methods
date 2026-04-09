# A. VFI components

# B. Backward logic 

# C. Forward logic  

# D. Solving 

# -----------------------------------------------
# A. VFI components 
# -----------------------------------------------

# 1. Static policies (MIT)
function fnStaticPoliciesMIT!(params, mit_endo, A⃗, λ⃗)
    
    # A. Unpacking business 
    @unpack α, δ, θ, z⃗, λ, a⃗,Tᴹᴵᵀ = params
    
    for it in 1:Tᴹᴵᵀ
        A1 = @. α * A⃗[it] * z⃗ / (mit_endo.rₜ[it] + δ)
        A2 = θ / α * (mit_endo.rₜ[it] + δ) / mit_endo.wₜ[it]
        A3 = A2^(θ)

        # Cleaned up manual dots
        kˣ                      = @. (A1 * A3)^(1 / (1 - α - θ))
        @. mit_endo.𝐤[:,:,it]   = min(kˣ, λ⃗[it] * a⃗')
        @. mit_endo.𝕀ᶜ[:,:,it]  = kˣ >= λ⃗[it] * a⃗'
        @. mit_endo.𝐥[:,:,it]   = (θ * A⃗[it] * z⃗ * mit_endo.𝐤[:,:,it]^α / mit_endo.wₜ[it])^(1 / (1 - θ))
        @. mit_endo.Π[:,:,it]   = A⃗[it] * z⃗ * mit_endo.𝐤[:,:,it]^α * mit_endo.𝐥[:,:,it]^θ - mit_endo.wₜ[it] * mit_endo.𝐥[:,:,it] - (mit_endo.rₜ[it] + δ) * mit_endo.𝐤[:,:,it]
    end 
end 