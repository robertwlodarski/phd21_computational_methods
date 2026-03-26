# Content:

# A. VFI components
# 1. Utility function 
# 2. Static policy functions 
# 3. Initial VFI guesses 
# 4. Solve for assets 
# 5. VFI

# B. Forward iteration 
# 1. Job desctruction
# 2. Update labour market (unemployment and matching)
# 3. Update the marginal distribution

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
    @unpack c̲, a⃗, β = params 

    # B. Cash variables
    Cash_w          = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ
    Cash_e          = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz,ia] - endo.τₜ

    # C. Bounds 
    Lower           = 0.0
    Upper_w         = Cash_w - c̲
    Upper_e         = Cash_e - c̲

    # D. Computation of assets: Worker 
    if Cash_w <= c̲
        A_w         = a⃗[1]
    else
        Obj(A_next) = -(fnUtility(Cash_w-A_next,params) + β * RHS_spline(A_next))
        res_w       = optimize(Obj,Lower,Upper_w,Brent())
        A_w         = Optim.minimizer(res_w)
    end 

    # E. Computation of assets: Entrepreneurs 
    if Cash_e <= c̲
        A_e         = a⃗[1]
    else
        Obj(A_next) = -(fnUtility(Cash_e-A_next,params) + β * RHS_spline(A_next))
        res_e       = optimize(Obj,Lower,Upper_e,Brent())
        A_e         = Optim.minimizer(res_e)
    end 
    return A_w, A_e
end

# 5. Value function iteration 
function fnVFI!(params, endo)

    # A. Unpacking business
     @unpack ψ, μ⃗, δᵛᶠⁱ, 𝒾̄ᵛᶠⁱ, z⃗, c̲, a⃗, β = params 

    # B. Static policies
    fnStaticPolicies!(params, endo)

    # C. Initial guess 
    fnInitialVFIGuess!(params, endo)

    # D. Prepare the VFI loop 
    εᵛᶠⁱ        = 1.0 
    𝓃ᵛᶠⁱ        = 1

    while (εᵛᶠⁱ > δᵛᶠⁱ && 𝓃ᵛᶠⁱ < 𝒾̄ᵛᶠⁱ)

        # D1. Update the expected VF and open the loop for productivity 
        endo.𝔼𝐕     .= ψ .* endo.𝐕 .+ (1 - ψ) .* (μ⃗' *  endo.𝐕)
        𝐕ᵖʳᵉᵛ       = copy(endo.𝐕)
        for iz in eachindex(z⃗)
        
            # B1. Prepare splines 
            ℑᶠ              = Spline1D(a⃗,(@views endo.𝔼𝐕[iz, :]); k=1, bc="extrapolate")
            
            # B2. Solve for each a and z 
            for ia in eachindex(a⃗)

                # I. Find assets
                aʷ, aᵉ          = fnFindAssets(iz, ia, ℑᶠ, params, endo)
                
                # II. VFs for different occupations
                cʷ              = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ - aʷ
                cᵉ              = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz,ia] - endo.τₜ - aᵉ
                endo.𝐕ᵂ[iz,ia]  = fnUtility(cʷ, params) + β * ℑᶠ(aʷ)
                endo.𝐕ᴱ[iz,ia]  = fnUtility(cᵉ, params) + β * ℑᶠ(aᵉ)

                # III. Occupational decision and updated VF
                endo.𝐨[iz,ia]   = (endo.𝐕ᴱ[iz,ia] >= endo.𝐕ᵂ[iz,ia])
                endo.𝐚[iz,ia]   = endo.𝐨[iz,ia] ? aᵉ : aʷ
                endo.𝐜[iz,ia]   = max(endo.𝐨[iz,ia] ? cᵉ : cʷ, c̲)
                endo.𝐕[iz,ia]   = max(endo.𝐕ᴱ[iz,ia], endo.𝐕ᵂ[iz,ia])
            end 
        end

        # D2. Error and update 
        εᵛᶠⁱ    = maximum(abs.(𝐕ᵖʳᵉᵛ.-endo.𝐕))
        𝓃ᵛᶠⁱ    += 1
    end 
end 

# -----------------------------------------------
# B. Forward iteration 
# -----------------------------------------------

# 1. Job desctruction
function fnJobDestruction!(params, endo)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, z⃗, a⃗, l⃗ = params 

    # B. Initialise 
    endo.S  = 0.0
    endo.D  = 0.0

    # C. Loop 
    for iz in eachindex(z⃗)
        for ia in eachindex(a⃗)
            for il in eachindex(l⃗)

                # C1. Compute key elements 
                Mass            = endo.G̃[iz, ia, il]
                LabDemand       = endo.𝐨[iz,ia] ? endo.𝐥[iz,ia] : 0.0
                JobsDestroyed   = max(0.0, l⃗[il] - LabDemand)
                Switchers       = (l⃗[il] > 0) && (endo.𝐨[iz,ia] == false)

                # C2. Start adding up 
                endo.S          += Mass * Switchers
                endo.D          += Mass * JobsDestroyed
            end 
        end 
    end 
    endo.JD     =   endo.S + endo.D
end

# 2. Update labour market (unemployment and matching)
function fnUpdateLabourMarket!(params, endo)

    # A. Unpacking business 
    @unpack γ = params 

    # B. Compute the key values
    endo.U      = sum(@views(endo.G[:,:,:,2]))
    endo.M      = γ * (endo.U + endo.JD)
end 

# 3. Update the marginal distribution
function fnUpdateMarginalDistribution(params, endo)

    # A. Unpacking business 
    @unpack z⃗, a⃗, l⃗, δᵈⁱˢᵗ = params 

    # B. Starting the loop 
    εᵈⁱˢᵗ       = 1.0
    while (εᵈⁱˢᵗ > δᵈⁱˢᵗ)
        for ia in eachindex(a⃗)
            for il in eachindex(l⃗)
                # I. The same productivity 
                # (a) Next asset mass 

                # (b) Labour mass 

                # (c) Update

                # II. With switching productivity 
                for iz in eachindex(z⃗)

                    # (a) Next asset mass 

                    # (b) Labour mass 

                    # (c) Update
                end 
            end 
        end 
    end 
end 

# 5. Forward iteration 
function fnDistributionIteration!(params, endo)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, δᵈⁱˢᵗ = params
    
    # B. Initial distribution guess 
    endo.G      .= ones(Nᶻ, Nᵃ, Nˡ,2) ./ (Nᶻ * Nᵃ * Nˡ * 2)
    endo.G̃      .= dropdims(sum(endo.G,dims=4),dims=4)
    
end 
