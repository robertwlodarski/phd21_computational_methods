# A. Backward logic 
# 1. Static policies
# 2. Last period 
# 3. Solve for assets 
# 4. Backward VF loop

# B. Forward simulation  
# 1. Job destruction
# 2. Update labour market
# 3. Compute flows
# 4. Distribution iteration
# 5. Aggregate states

# C. Solving 

# -----------------------------------------------
# A. Backward logic 
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

# 2. Last period (MIT)
function fnLastPeriodMIT!(params, mit_endo, ss_endo)

    # A. Unpacking business 
    @unpack Tᴹᴵᵀ = params 

    # B. Value functions 
    mit_endo.𝐕[:,:,Tᴹᴵᵀ]    .= ss_mit_endo.𝐕
    mit_endo.𝔼𝐕[:,:,Tᴹᴵᵀ]   .= ss_mit_endo.𝔼𝐕
    mit_endo.𝐕ᵂ[:,:,Tᴹᴵᵀ]   .= ss_mit_endo.𝐕ᵂ
    mit_endo.𝐕ᴱ[:,:,Tᴹᴵᵀ]   .= ss_mit_endo.𝐕ᴱ

    # C. Policy functions 
    mit_endo.𝐨[:,:,Tᴹᴵᵀ]    .= ss_mit_endo.𝐨
    mit_endo.𝐤[:,:,Tᴹᴵᵀ]    .= ss_mit_endo.𝐤
    mit_endo.𝐚[:,:,Tᴹᴵᵀ]    .= ss_mit_endo.𝐚
    mit_endo.𝐜[:,:,Tᴹᴵᵀ]    .= ss_mit_endo.𝐜
    mit_endo.𝐥[:,:,Tᴹᴵᵀ]    .= ss_mit_endo.𝐥
end 

# 3. Solve for assets 
# A. A struct to hold the explicitly typed variables
struct AssetObjective{S, P}
    cash::Float64
    spline::S
    params::P
end

# B. Make the struct callable (this replaces the anonymous function)
function (obj::AssetObjective)(ap)
    return -(fnUtility(obj.cash - ap, obj.params) + obj.spline(ap))
end

# C. Start the function 
function fnFindAssetsMIT(iz, ia, it, exp_spline, params, mit_endo)

    # D. Unpacking parameters
    @unpack c̲, a⃗, β = params

    # E. Helper elements 
    Cash_w = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.wₜ[it] - mit_endo.τₜ[it]
    Cash_e = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.Π[iz, ia,it] - mit_endo.τₜ[it]
    Lower   = a⃗[1]
    Upper_w = min(Cash_w - c̲, a⃗[end])
    Upper_e = min(Cash_e - c̲, a⃗[end])

    # F. Computation of assets: Worker 
    if Upper_w <= Lower
        A_w = a⃗[1]
        𝐕ᵂ  = fnUtility(Cash_w - a⃗[1], params) + exp_spline(a⃗[1])
        𝔼𝐕ᵂ = exp_spline(A_w)
    else
        obj_w = AssetObjective(Cash_w, exp_spline, params)
        ℜʷ    = optimize(obj_w, Lower, Upper_w)
        A_w   = Optim.minimizer(ℜʷ)
        𝐕ᵂ    = -Optim.minimum(ℜʷ)
        𝔼𝐕ᵂ   = exp_spline(A_w)
    end

    # G. Computation of assets: Entrepreneurs 
    if Upper_e <= Lower
        A_e = a⃗[1]
        𝐕ᴱ  = fnUtility(Cash_e - a⃗[1], params) + exp_spline(a⃗[1])
        𝔼𝐕ᴱ = exp_spline(A_e)
    else
        obj_e = AssetObjective(Cash_e, exp_spline, params)
        ℜᴱ    = optimize(obj_e, Lower, Upper_e)
        A_e   = Optim.minimizer(ℜᴱ)
        𝐕ᴱ    = -Optim.minimum(ℜᴱ)
        𝔼𝐕ᴱ   = exp_spline(A_e)
    end
    
    # H. Returning business  
    return A_w, A_e, 𝐕ᵂ, 𝐕ᴱ, 𝔼𝐕ᵂ, 𝔼𝐕ᴱ
end

# 4. Backward VF loop
function fnBackwardInductionMIT!(params, mit_endo, ss_endo, A⃗, λ⃗)

    # A. Unpacking business 
    @unpack Tᴹᴵᵀ,a⃗,z⃗,μ⃗,β,ψ,c̲,Nᶻ,Nᵃ = params 

    # B. Prepare the setting 
    fnStaticPoliciesMIT!(params, mit_endo, A⃗, λ⃗)

    # C. Impose the steady state in the last period
    fnLastPeriodMIT!(params, mit_endo, ss_endo)

    # D. Start the loop 
    𝔼𝐕ᵂ     = zeros(Nᶻ, Nᵃ,Tᴹᴵᵀ) 
    𝔼𝐕ᴱ     = zeros(Nᶻ, Nᵃ,Tᴹᴵᵀ)

    for it in Tᴹᴵᵀ-1:(-1):1

        # D1. Spline business 
        ℑᴱ                      = [Spline1D(a⃗, mit_endo.𝐕ᴱ[iz, :,it+1]; k=1, bc="nearest") for iz in eachindex(z⃗)]
        ℑᵂ                      = [Spline1D(a⃗, mit_endo.𝐕ᵂ[iz, :,it+1]; k=1, bc="nearest") for iz in eachindex(z⃗)]
        𝔼_max(a)                = sum(μ⃗[jz] * max(ℑᴱ[jz](a), ℑᵂ[jz](a)) for jz in eachindex(z⃗))
        ℑ(iz, a)                = β * (ψ * max(ℑᴱ[iz](a), ℑᵂ[iz](a)) + (1 - ψ) * 𝔼_max(a))

        # D2. Open the inner loops 
        Threads.@threads for ia in eachindex(a⃗)
            @inbounds for iz in eachindex(z⃗)

                # I. Find assets
                ℑᶻ(a)   = ℑ(iz, a)
                mit_endo.𝐚ʷ[iz, ia, it], mit_endo.𝐚ᵉ[iz, ia, it], mit_endo.𝐕ᵂ[iz, ia, it], mit_endo.𝐕ᴱ[iz, ia, it], 𝔼𝐕ᵂ[iz, ia, it], 𝔼𝐕ᴱ[iz, ia, it] = fnFindAssetsMIT(iz, ia, it, ℑᶻ, params, mit_endo)
                
                # II. VFs for different occupations
                mit_endo.𝐜ʷ[iz, ia, it] = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.wₜ[it] - mit_endo.τₜ[it] - mit_endo.𝐚ʷ[iz, ia, it]
                mit_endo.𝐜ᵉ[iz, ia, it] = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.Π[iz, ia, it] - mit_endo.τₜ[it] - mit_endo.𝐚ᵉ[iz, ia, it]

                # III. Occupational decision and updated VF
                mit_endo.𝐨[iz, ia, it]  = (mit_endo.𝐕ᴱ[iz, ia, it] >= mit_endo.𝐕ᵂ[iz, ia, it])
                mit_endo.𝐚[iz, ia, it]  = mit_endo.𝐨[iz, ia, it] ? mit_endo.𝐚ᵉ[iz, ia, it] : mit_endo.𝐚ʷ[iz, ia, it]
                mit_endo.𝐜[iz, ia, it]  = max(mit_endo.𝐨[iz, ia, it] ? mit_endo.𝐜ᵉ[iz, ia, it] : mit_endo.𝐜ʷ[iz, ia, it], c̲)
                mit_endo.𝐕[iz, ia, it]  = max(mit_endo.𝐕ᴱ[iz, ia, it], mit_endo.𝐕ᵂ[iz, ia, it])
            end 
        end 
    end 
end 

# -----------------------------------------------
# B. Forward simulation 
# -----------------------------------------------

# 1. Job destruction (MIT)
function fnJobDestructionMIT!(params, mit_endo, it)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, z⃗, a⃗, l⃗ = params

    # B. Initialise 
    mit_endo.S[it]  = 0.0
    mit_endo.D[it]  = 0.0

    # C. Loop 
    for iz in eachindex(z⃗)
        for ia in eachindex(a⃗)
            for il in eachindex(l⃗)

                # C1. Compute key elements 
                Mass            = mit_endo.g[iz, ia, il, 1, it]
                LabDemand       = mit_endo.𝐨[iz, ia, it] ? mit_endo.𝐥[iz, ia, it] : 0.0
                JobsDestroyed   = max(0.0, l⃗[il] - LabDemand)
                Switchers       = (l⃗[il] > 0) && (mit_endo.𝐨[iz, ia, it] == false)

                # C2. Start adding up 
                mit_endo.S[it] += Mass * Switchers
                mit_endo.D[it] += Mass * JobsDestroyed
            end
        end
    end
    mit_endo.JD[it] = mit_endo.S[it] + mit_endo.D[it]
end

# 2. Update labour market (MIT)
function fnUpdateLabourMarketMIT!(params, mit_endo, it)

    # A. Unpacking business 
    @unpack γ = params

    # B. Compute the key values
    mit_endo.U[it]  = sum(@views(mit_endo.g[:, :, :, 2, it]))
    mit_endo.M[it]  = γ * (mit_endo.U[it] + mit_endo.JD[it])
    mit_endo.W[it]  = sum(@views(mit_endo.g[:, :, 1, 1, it]))
    mit_endo.E[it]  = 1.0 - mit_endo.U[it] - mit_endo.W[it]
end

# 3. Compute flows (MIT)
function fnComputeFlowsMIT(iz, ia, il, iu, it, params, mit_endo)

    # A. Unpacking business
    @unpack l⃗ = params
    # Legend: 
    # iu    = 1 -> employed, entrepreneur 
    # iu    = 2 -> unemployed 
    # E     = Entrepreneurs 
    # W     = Employed workers 
    # U     = Unemployed workers

    # B. Flow indicators 
    JFR = min(1.0, max(mit_endo.M[it] / (mit_endo.U[it] + mit_endo.JD[it]), 0.0))
    JDR = min(1.0, max(0.0, mit_endo.D[it] / mit_endo.W[it]))
    WU  = 0.0 + (l⃗[il] == l⃗[1]) * JDR * (1 - JFR) * (iu == 1) * (mit_endo.𝐨[iz, ia, it] == false)
    WE  = 0.0 + (l⃗[il] == l⃗[1]) * (mit_endo.𝐨[iz, ia, it] == true) * (iu == 1)
    EU  = 0.0 + (l⃗[il] > l⃗[1]) * (mit_endo.𝐨[iz, ia, it] == false) * (iu == 1) * (1 - JFR)
    EW  = 0.0 + (l⃗[il] > l⃗[1]) * (mit_endo.𝐨[iz, ia, it] == false) * (iu == 1) * JFR
    UE  = 0.0 + (mit_endo.𝐨[iz, ia, it] == true) * (iu == 2)
    UW  = 0.0 + (mit_endo.𝐨[iz, ia, it] == false) * (iu == 2) * JFR
    EE  = 0.0 + (l⃗[il] > l⃗[1]) * (mit_endo.𝐨[iz, ia, it] == true) * (iu == 1)
    WW  = 0.0 + (l⃗[il] == l⃗[1]) * (mit_endo.𝐨[iz, ia, it] == false) * (iu == 1) * ((1 - JDR) + (JDR * JFR))
    UU  = 0.0 + (mit_endo.𝐨[iz, ia, it] == false) * (iu == 2) * (1 - JFR)

    # C. Combine indicators to account for flows that matter
    F1t2 = WU + EU           # To unemployment
    F2t1 = UW + UE           # From unemployment 
    F2t2 = UU
    F1t1 = EE + WW + WE + EW
    return F1t2, F2t1, F2t2, F1t1
end

# 4. Distribution iteration (MIT)
function fnDistributionIterationMIT!(params, mit_endo, ss_endo)

    # A. Unpacking business 
    @unpack z⃗, a⃗, l⃗, ψ, μ⃗, Nᶻ, Nᵃ, Nˡ, Nᵘ, Tᴹᴵᵀ = params

    # B. Starting the loop for a PDF at different times 
    B                           = zeros(Nᵃ, Nˡ, Nᵘ,Tᴹᴵᵀ)
    gⁿᵉˣᵗ                       = zeros(Nᶻ, Nᵃ, Nˡ, Nᵘ,Tᴹᴵᵀ)
    @. mit_endo.g[:,:,:,:, 1]   = @views(ss_endo.g)
    gᵖʳᵉᵛ                       = zeros(Nᶻ, Nᵃ, Nˡ, Nᵘ,Tᴹᴵᵀ)
    @. gᵖʳᵉᵛ[:,:,:,:, 1]        = @views(ss_endo.g)

    # C. Precompute the invariant elements 
    # I. Initialisation 
    ibᵃ     = zeros(Int, Nᶻ, Nᵃ, Tᴹᴵᵀ)
    iuᵃ     = zeros(Int, Nᶻ, Nᵃ, Tᴹᴵᵀ)
    wᵃ      = zeros(Float64, Nᶻ, Nᵃ, Tᴹᴵᵀ)
    ibˡ     = zeros(Int, Nᶻ, Nᵃ, Tᴹᴵᵀ)
    iuˡ     = zeros(Int, Nᶻ, Nᵃ, Tᴹᴵᵀ)
    wˡ      = zeros(Float64, Nᶻ, Nᵃ, Tᴹᴵᵀ)

    # II. Loop 
    @inbounds for it in 1:Tᴹᴵᵀ
        for ia in eachindex(a⃗)
            for iz in eachindex(z⃗)

                # (a) Next asset mass 
                Aₜ              = mit_endo.𝐚[iz, ia, it]
                ibᵃ[iz, ia, it] = clamp(searchsortedlast(a⃗, Aₜ), 1, length(a⃗) - 1)
                iuᵃ[iz, ia, it] = ibᵃ[iz, ia, it] + 1
                wᵃ[iz, ia, it]  = clamp((a⃗[iuᵃ[iz, ia, it]] - Aₜ) / (a⃗[iuᵃ[iz, ia, it]] - a⃗[ibᵃ[iz, ia, it]]), 0.0, 1.0)
                
                # (b) Labour mass 
                Lₜ              = mit_endo.𝐨[iz, ia, it] ? mit_endo.𝐥[iz, ia, it] : 0.0
                ibˡ[iz, ia, it] = clamp(searchsortedlast(l⃗, Lₜ), 1, length(l⃗) - 1)
                iuˡ[iz, ia, it] = ibˡ[iz, ia, it] + 1
                wˡ[iz, ia, it]  = clamp((l⃗[iuˡ[iz, ia, it]] - Lₜ) / (l⃗[iuˡ[iz, ia, it]] - l⃗[ibˡ[iz, ia, it]]), 0.0, 1.0)
            end
        end
    end
    
    # D. Start the largest loop 
    @inbounds for it in 1:Tᴹᴵᵀ-1

        # I. Update the endogenous functions for time t
        fnJobDestructionMIT!(params, mit_endo, it)
        fnUpdateLabourMarketMIT!(params, mit_endo, it)

        for iu in 1:Nᵘ
            for il in eachindex(l⃗)
                for ia in eachindex(a⃗)
                    for iz in eachindex(z⃗)

                        # II. Common terms
                        Mass = gᵖʳᵉᵛ[iz, ia, il, iu, it]
                        if Mass == 0.0
                            continue
                        end

                        # III. Update the flows 
                        # (c) Get flows that matter
                        f¹², f²¹, f²², f¹¹ = fnComputeFlowsMIT(iz, ia, il, iu, it, params, mit_endo)
                        FlowU = f¹² + f²²
                        FlowO = f²¹ + f¹¹

                        # IV. With switching productivity
                        # To unemployment (iu=2)
                        B[ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 2, it]      += Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowU
                        B[ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 2, it]      += Mass * wᵃ[iz,ia, it] * (1 - wˡ[iz,ia,it]) * FlowU
                        B[iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 2, it]      += Mass * (1 - wᵃ[iz,ia,it]) * wˡ[iz,ia,it] * FlowU
                        B[iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 2, it]      += Mass * (1 - wᵃ[iz,ia,it]) * (1 - wˡ[iz,ia,it]) * FlowU
                        # To employment/entrepreneurship (iu=1)
                        B[ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 1, it]      += Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowO
                        B[ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 1, it]      += Mass * wᵃ[iz,ia,it] * (1 - wˡ[iz,ia,it]) * FlowO
                        B[iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 1, it]      += Mass * (1 - wᵃ[iz,ia,it]) * wˡ[iz,ia,it] * FlowO
                        B[iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 1, it]      += Mass * (1 - wᵃ[iz,ia,it]) * (1 - wˡ[iz,ia,it]) * FlowO

                        # V. The same productivity
                        # (c) Update (asset, labour) = (BB, BU, UB, UU)
                        # To unemployment 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 2, it]      += ψ * Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowU
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 2, it]      += ψ * Mass * wᵃ[iz,ia,it] * (1 - wˡ[iz,ia,it]) * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 2, it]      += ψ * Mass * (1 - wᵃ[iz,ia, it]) * wˡ[iz,ia, it] * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 2, it]      += ψ * Mass * (1 - wᵃ[iz,ia, it]) * (1 - wˡ[iz,ia,it]) * FlowU
                        # To others 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 1, it]      += ψ * Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowO
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 1, it]      += ψ * Mass * wᵃ[iz,ia,it] * (1 - wˡ[iz,ia,it]) * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 1, it]      += ψ * Mass * (1 - wᵃ[iz,ia,it]) * wˡ[iz,ia,it] * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 1, it]      += ψ * Mass * (1 - wᵃ[iz,ia,it]) * (1 - wˡ[iz,ia,it]) * FlowO
                    end
                end
            end
        end

        # VI. Add the mass to the updated productivity states
        @inbounds for iu in 1:2, il in 1:Nˡ ,ia in 1:Nᵃ, iz in 1:Nᶻ
            gⁿᵉˣᵗ[iz, ia, il, iu, it] += (1 - ψ) * μ⃗[iz] * B[ia, il, iu, it]
        end

        # VII. Update 
        @. mit_endo.g[:,:,:,:,it+1] = @views(gⁿᵉˣᵗ[:,:,:,:,it])
        @. gᵖʳᵉᵛ[:,:,:,:,it+1]      = @views(mit_endo.g[:,:,:,:,it+1])
    end
end

# 5. Aggregate states 
function fnAggregateStatesMIT!(params, mit_endo, ss_endo)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, a⃗, Tᴹᴵᵀ = params

    # B. Iterate to find the distribution and its marginal 
    fnDistributionIterationMIT!(params, mit_endo, ss_endo)

    for it in 1:Tᴹᴵᵀ
        # C. Weight for the marginal z × a distribution
        ωᵃ              = dropdims(sum(@views(mit_endo.g[:,:,:,:,it]), dims=(3, 4)), dims=(3, 4))

        # D. Capital demand and supply 
        mit_endo.Kˢ[it] = sum(ωᵃ .* a⃗')
        mit_endo.Kᵈ[it] = sum(ωᵃ .* @views(mit_endo.𝐤[:,:,it]) .*  @views(mit_endo.𝐨[:,:,it]))

        # E. Labour demand and supply 
        mit_endo.Lˢ[it] = mit_endo.W[it]
        mit_endo.Lᵈ[it] = sum(ωᵃ .* @views(mit_endo.𝐥[:,:,it]) .* @views(mit_endo.𝐨[:,:,it]))
    end 
end