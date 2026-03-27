# Content:

# A. VFI components
# 1. Utility function 
# 2. Static policy functions 
# 3. Initial VFI guesses 
# 4. Solve for assets 
# 5. VFI

# B. Aggregating 
# 1. Job destruction
# 2. Update labour market (unemployment and matching)
# 3. Compute flows 
# 4. Forward iteration
# 5. Aggregate states 

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
# B. Aggregating 
# -----------------------------------------------

# 1. Job destruction
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
                Mass            = endo.g[iz, ia, il, 1]
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
    endo.U      = sum(@views(endo.g[:,:,:,2]))
    endo.M      = γ * (endo.U + endo.JD)
    endo.W      = sum(@views(endo.g[:,:,1,1]))
    endo.E      = 1.0 - endo.U - endo.W
end 

# 3. Compute flows 
function fnComputeFlows(iz,ia,il,iu, params, endo)

    # A. Unpacking business
    @unpack l⃗ = params 
    # Legend: 
    # iu    = 1 -> employed, entrepreneur 
    # iu    = 2 -> unemployed 
    # E     = Entrepreneurs 
    # W     = Employed workers 
    # U     = Unemployed workers

    # B. Flow indicators 
    JFR = min(1.0, max(endo.M / (endo.U + endo.JD), 0.0))
    JDR = min(1.0, max(0.0, endo.JD / endo.W))
    WU  = 0.0 + (l⃗[il] == l⃗[1]) * JDR * (1 - JFR) * (iu == 1) * (endo.𝐨[iz, ia] == false)
    WE  = 0.0 + (l⃗[il] == l⃗[1]) * (endo.𝐨[iz, ia] == true) * (iu == 1)
    EU  = 0.0 + (l⃗[il] > l⃗[1]) * (endo.𝐨[iz, ia] == false) * (iu == 1) * (1 - JFR)
    EW  = 0.0 + (l⃗[il] > l⃗[1]) * (endo.𝐨[iz, ia] == false) * (iu == 1) * JFR
    UE  = 0.0 + (endo.𝐨[iz, ia] == true) * (iu == 2)
    UW  = 0.0 + (endo.𝐨[iz, ia] == false) * (iu == 2) * JFR
    EE  = 0.0 + (l⃗[il] > l⃗[1]) * (endo.𝐨[iz, ia] == true) * (iu == 1)
    WW  = 0.0 + (endo.𝐨[iz, ia] == false) * (iu == 1) * ((1 - JDR) + (JDR * JFR))
    UU  = 0.0 + (endo.𝐨[iz, ia] == false) * (iu == 2) * (1 - JFR)

    # C. Combine indicators to account for flows that matter
    F1t2 = WU + EU           # To unemployment
    F2t1 = UW + UE           # From unemployment 
    F2t2 = UU
    F1t1 = EE + WW + WE + EW 
    return F1t2, F2t1, F2t2, F1t1
end 

# 4. Forward iteration 
function fnForwardIteration!(params, endo)

    # A. Unpacking business 
    @unpack z⃗, a⃗, l⃗, δᵈⁱˢᵗ, ψ,μ⃗ = params 

    # B. Starting the loop for a PDF
    εᵈⁱˢᵗ       = 1.0
    gⁿᵉˣᵗ       = zeros(Nᶻ, Nᵃ, Nˡ,2)     
    while (εᵈⁱˢᵗ > δᵈⁱˢᵗ)
        
        # C. Aggregate states 
        fnJobDestruction!(params, endo)         # Getting JD, S 
        fnUpdateLabourMarket!(params, endo)     # Getting U, M, W, E  

        for ia in eachindex(a⃗)
            for il in eachindex(l⃗)
                for iz in eachindex(z⃗)
                    for iu in 1:2
                        # I. Get the key flows 
                       

                        # II. Common terms
                        Mass        = endo.g[iz,ia,il,iu] 
                        if Mass == 0.0
                            continue 
                        end
                        # (a) Next asset mass 
                        Aₜ                  = endo.𝐚[iz,ia]
                        ibᵃ                 = clamp(searchsortedlast(a⃗,Aₜ),1,length(a⃗)-1)
                        iuᵃ                 = ibᵃ + 1
                        wᵃ                  = clamp((a⃗[iuᵃ]-Aₜ) / (a⃗[iuᵃ]-a⃗[ibᵃ]),0.0,1.0)
                        # (b) Labour mass 
                        Lₜ                  = endo.𝐨[iz,ia] ? endo.𝐥[iz,ia] : 0.0
                        ibˡ                 = clamp(searchsortedlast(l⃗,Lₜ),1,length(l⃗)-1)
                        iuˡ                 = ibˡ + 1
                        wˡ                  = clamp((l⃗[iuˡ]-Lₜ) / (l⃗[iuˡ]-l⃗[ibˡ]),0.0,1.0)
                        # (c) Get flows that matter
                            f¹², f²¹ , f²², f¹¹     = fnComputeFlows(iz,ia,il,iu, params, endo)
                            FlowU                   = f¹² + f²²
                            FlowO                   = f²¹ + f¹¹

                        # III. With switching productivity
                        for izp in eachindex(z⃗)
                            # (d) Update (asset, labour) = (BB, BU, UB, UU)
                            # To unemployment 
                            gⁿᵉˣᵗ[izp,ibᵃ,ibˡ,2]   += (1 - ψ) * μ⃗[izp] * Mass * wᵃ * wˡ * FlowU
                            gⁿᵉˣᵗ[izp,ibᵃ,iuˡ,2]   += (1 - ψ) * μ⃗[izp] * Mass * wᵃ * (1 - wˡ) * FlowU
                            gⁿᵉˣᵗ[izp,iuᵃ,ibˡ,2]   += (1 - ψ) * μ⃗[izp] * Mass * (1- wᵃ) * wˡ * FlowU
                            gⁿᵉˣᵗ[izp,iuᵃ,iuˡ,2]   += (1 - ψ) * μ⃗[izp] * Mass * (1- wᵃ) * (1 - wˡ) * FlowU
                            # To others 
                            gⁿᵉˣᵗ[izp,ibᵃ,ibˡ,1]   += (1 - ψ) * μ⃗[izp] * Mass * wᵃ * wˡ * FlowO
                            gⁿᵉˣᵗ[izp,ibᵃ,iuˡ,1]   += (1 - ψ) * μ⃗[izp] * Mass * wᵃ * (1 - wˡ) * FlowO
                            gⁿᵉˣᵗ[izp,iuᵃ,ibˡ,1]   += (1 - ψ) * μ⃗[izp] * Mass * (1- wᵃ) * wˡ * FlowO
                            gⁿᵉˣᵗ[izp,iuᵃ,iuˡ,1]   += (1 - ψ) * μ⃗[izp] * Mass * (1- wᵃ) * (1 - wˡ) * FlowO
                        end 

                        # IV. The same productivity
                        # (c) Update (asset, labour) = (BB, BU, UB, UU)
                        # To unemployment 
                        gⁿᵉˣᵗ[iz, ibᵃ, ibˡ, 2] += ψ * Mass * wᵃ * wˡ * FlowU
                        gⁿᵉˣᵗ[iz, ibᵃ, iuˡ, 2] += ψ * Mass * wᵃ * (1 - wˡ) * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ, ibˡ, 2] += ψ * Mass * (1 - wᵃ) * wˡ * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ, iuˡ, 2] += ψ * Mass * (1 - wᵃ) * (1 - wˡ) * FlowU
                        # To others 
                        gⁿᵉˣᵗ[iz, ibᵃ, ibˡ, 1] += ψ * Mass * wᵃ * wˡ * FlowO
                        gⁿᵉˣᵗ[iz, ibᵃ, iuˡ, 1] += ψ * Mass * wᵃ * (1 - wˡ) * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ, ibˡ, 1] += ψ * Mass * (1 - wᵃ) * wˡ * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ, iuˡ, 1] += ψ * Mass * (1 - wᵃ) * (1 - wˡ) * FlowO
                    end 
                end 
            end 
        end 

        # VI. Updating business 
        εᵈⁱˢᵗ       = maximum(abs.(gⁿᵉˣᵗ .- endo.g))
        endo.g      .= gⁿᵉˣᵗ
        fill!(gⁿᵉˣᵗ, 0.0)
    end 
end 

# 5. Aggregate states 
function fnAggregateStates!(params, endo)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, δᵈⁱˢᵗ = params
    
    # B. Initial PDF guess
    endo.g      .= ones(Nᶻ, Nᵃ, Nˡ,2) ./ (Nᶻ * Nᵃ * Nˡ * 2)

    # C. Iterate to find the distribution and its marginal 
    fnForwardIteration!(params, endo)
    endo.g̃      .= dropdims(sum(endo.g,dims=4),dims=4)

    # D. Update the labour market 
    fnJobDestruction!(params, endo)
    fnUpdateLabourMarket!(params, endo)

    # E. Weight for the marginal z × a distribution
    ωᵃ      = dropdims(sum(endo.g, dims=(3,4)), dims=(3,4))

    # F. Capital demand and supply 
    endo.Kˢ = sum(ωᵃ .* endo.𝐚)     
    endo.Kᵈ = sum(ωᵃ .* endo.𝐤 .* endo.𝐨) 

    # G. Labour demand and supply 
    endo.Lˢ = endo.W 
    endo.Lᵈ = sum(ωᵃ .* endo.𝐥 .* endo.𝐨)
end 

# -----------------------------------------------
# 3. Solving  
# -----------------------------------------------