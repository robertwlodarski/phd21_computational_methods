# Content:
# A. VFI components
# 1. Utility function 
# 2. Static policy functions 
# 3. Initial VFI guesses 

# B. Forward iteration 
# 1. 


# -----------------------------------------------
# A. VFI components 
# -----------------------------------------------

# 1. Utility function 
function fnUtility(C,params)
    
    # A. Unpacking business 
    @unpack σ   = params 

    # B. Compute the utility 
    return C.^(1-σ) ./ (1-σ) 
end 

# 2. Static policy functions  
function fnStaticPolicies!(params, endo)

    # A. Unpacking business 
    @unpack α, δ, θ, A, z⃗, λ, a⃗ = params 

    # B. Simplifying notation 
    A1      = α * A * z⃗ ./ (endo.rₜ + δ)
    A2      = θ / α * (endo.rₜ + δ) / endo.wₜ
    A3      = A2^(θ)

    # C. Getting the capital 
    kˣ      = (A1 .* A3).^(1 / (1 - α - θ))
    endo.𝐤  .= min.(kˣ,λ .* a⃗')
    endo.𝕀ᶜ .= (kˣ .>= λ .* a⃗')

    # D. Getting the labour 
    endo.𝐥  .= (θ .* A .* z⃗ .* endo.𝐤.^α ./ endo.wₜ).^(1 / (1 - θ))

    # E. Compute the associated profit 
    endo.Π .= A .* z⃗ .* (endo.𝐤).^α .* (endo.𝐥).^θ .- endo.wₜ .* endo.𝐥 .- (endo.rₜ + δ) .* endo.𝐤
end 

# 3. Initial VFI guesses 
function fnInitialVFIGuess!(params, endo)

    # A. Unpacking business 
    @unpack β, a⃗ = params 

    # B. Workers VF 
    endo.𝐕ᵂ     .= (1-β)^(-1) * fnUtility(endo.rₜ .* a⃗' .- endo.τₜ + endo.wₜ ,params)

    # C. Entrepreneurs VF 
    endo.𝐕ᴱ     .= (1-β)^(-1) * fnUtility(endo.rₜ .* a⃗' .- endo.τₜ + endo.Π ,params)

    # D. Overall VF 
    endo.𝐕      .= max.(endo.𝐕ᵂ,endo.𝐕ᴱ) 
    endo.𝐨      .= (endo.𝐕ᵂ .<= endo.𝐕ᴱ) 
end 

# 4. Solve for assets 
function fnFindAssets(iz, ia, RHS_spline, params, endo)

    # A. Unpack parameters 
   

    # B. Entrepreneurs assets 
    # Use the public econ code 

    return A_w, A_nw
end

# 5. Value function iteration 
function fnVFI!(params, endo)

    # A. Unpacking business
     @unpack ψ, μ⃗, δʳᵉᶠ, 𝒾̄ᵛᶠⁱ, z⃗ = params 

    # B. Static policies
    fnStaticPolicies!(params, endo)

    # C. Initial guess 
    fnInitialVFIGuess!(params, endo)

    # D. Prepare the VFI loop 
    εᵛᶠⁱ        = 1.0 
    𝓃ᵛᶠⁱ        = 1
    𝐕ⁿᵉˣᵗ       = copy(endo.𝐕)

    while (εᵛᶠⁱ > δʳᵉᶠ & 𝓃ᵛᶠⁱ < 𝒾̄ᵛᶠⁱ)

        # D1. Update the expected VF and open the loop for productivity 
        endo.𝔼𝐕     .= ψ .* endo.𝐕 .+ (1 - ψ) .* (μ⃗' *  endo.𝐕)
        for iz in eachindex(z⃗)
        
            # B1. Prepare splines 
            ℑᶠ          = Spline1D(a⃗,(@views endo.𝔼𝐕[iz, :]); k=1, bc="extrapolate")
            
            # B2. Solve for each a and z 
            for ia in eachindex(a⃗)

                # I. Find assets
                
                # II. VFs 

                # III. Occupational decision and updated VF
                
            end 

        end

        # D2. Error and update 

    end 

end 