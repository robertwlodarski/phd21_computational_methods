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
# 1. Convex updating for MIT shock
# 2. Labour market residual
# 3. Government balance 

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
    mit_endo.𝐕[:,:,Tᴹᴵᵀ]    .= ss_endo.𝐕
    mit_endo.𝔼𝐕[:,:,Tᴹᴵᵀ]   .= ss_endo.𝔼𝐕
    mit_endo.𝐕ᵂ[:,:,Tᴹᴵᵀ]   .= ss_endo.𝐕ᵂ
    mit_endo.𝐕ᴱ[:,:,Tᴹᴵᵀ]   .= ss_endo.𝐕ᴱ

    # C. Policy functions 
    mit_endo.𝐨[:,:,Tᴹᴵᵀ]    .= ss_endo.𝐨
    mit_endo.𝐤[:,:,Tᴹᴵᵀ]    .= ss_endo.𝐤
    mit_endo.𝐚[:,:,Tᴹᴵᵀ]    .= ss_endo.𝐚
    mit_endo.𝐜[:,:,Tᴹᴵᵀ]    .= ss_endo.𝐜
    mit_endo.𝐥[:,:,Tᴹᴵᵀ]    .= ss_endo.𝐥
    mit_endo.Π[:,:,Tᴹᴵᵀ]    .= ss_endo.Π
end 

# 3A. Solve for assets 
# A. A struct to hold the explicitly typed variables
# struct AssetObjective{S, P}
#     cash::Float64
#     spline::S
#     params::P
# end

# # B. Make the struct callable (this replaces the anonymous function)
# function (obj::AssetObjective)(ap)
#     return -(fnUtility(obj.cash - ap, obj.params) + obj.spline(ap))
# end

# # C. Start the function 
# function fnFindAssetsMIT(iz, ia, it, exp_spline, params, mit_endo)

#     # D. Unpacking parameters
#     @unpack c̲, a⃗, β = params

#     # E. Helper elements 
#     Cash_w = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.wₜ[it] - mit_endo.τₜ[it]
#     Cash_e = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.Π[iz, ia,it] - mit_endo.τₜ[it]
#     Lower   = a⃗[1]
#     Upper_w = min(Cash_w - c̲, a⃗[end])
#     Upper_e = min(Cash_e - c̲, a⃗[end])

#     # F. Computation of assets: Worker 
#     if Upper_w <= Lower
#         A_w = a⃗[1]
#         𝐕ᵂ  = fnUtility(Cash_w - a⃗[1], params) + exp_spline(a⃗[1])
#         𝔼𝐕ᵂ = exp_spline(A_w)
#     else
#         obj_w = AssetObjective(Cash_w, exp_spline, params)
#         ℜʷ    = optimize(obj_w, Lower, Upper_w)
#         A_w   = Optim.minimizer(ℜʷ)
#         𝐕ᵂ    = -Optim.minimum(ℜʷ)
#         𝔼𝐕ᵂ   = exp_spline(A_w)
#     end

#     # G. Computation of assets: Entrepreneurs 
#     if Upper_e <= Lower
#         A_e = a⃗[1]
#         𝐕ᴱ  = fnUtility(Cash_e - a⃗[1], params) + exp_spline(a⃗[1])
#         𝔼𝐕ᴱ = exp_spline(A_e)
#     else
#         obj_e = AssetObjective(Cash_e, exp_spline, params)
#         ℜᴱ    = optimize(obj_e, Lower, Upper_e)
#         A_e   = Optim.minimizer(ℜᴱ)
#         𝐕ᴱ    = -Optim.minimum(ℜᴱ)
#         𝔼𝐕ᴱ   = exp_spline(A_e)
#     end
    
#     # H. Returning business  
#     return A_w, A_e, 𝐕ᵂ, 𝐕ᴱ, 𝔼𝐕ᵂ, 𝔼𝐕ᴱ
# end

# 3B. Optimised grid search 
 function fnPureGridSearchMIT(Budget, params, a⃗, 𝔼V_array)
    
    # A. Unpacking business
    @unpack c̲, β = params
    
    # B. Initial values 
    best_val    = -Inf
    best_a      = a⃗[1]
    
    # C. Start the loop 
    @inbounds for i in eachindex(a⃗)
        
        # D. Savings and consumption 
        a_prime         = a⃗[i]
        c               = Budget - a_prime
        
        # D. Stop searching if consumption hits the minimum bound
        if c <= c̲
            v_bound             = fnUtility(c̲, params) + β * 𝔼V_array[i]
            if v_bound > best_val
                best_val        = v_bound
                best_a          = a_prime
            end
            break 
        end
        
        # E. Pure array lookup 
        v                       = fnUtility(c, params) + β * 𝔼V_array[i]
        if v > best_val
            best_val            = v
            best_a              = a_prime
        end
    end
    
    # F. Returning business 
    return best_a, best_val
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
    # 𝔼𝐕ᵂ     = zeros(Nᶻ, Nᵃ,Tᴹᴵᵀ) 
    # 𝔼𝐕ᴱ     = zeros(Nᶻ, Nᵃ,Tᴹᴵᵀ)
    
    for it in Tᴹᴵᵀ-1:(-1):1

        # D1. Spline business 
        # ℑᴱ                      = [Spline1D(a⃗, mit_endo.𝐕ᴱ[iz, :,it+1]; k=1, bc="nearest") for iz in eachindex(z⃗)]
        # ℑᵂ                      = [Spline1D(a⃗, mit_endo.𝐕ᵂ[iz, :,it+1]; k=1, bc="nearest") for iz in eachindex(z⃗)]
        # 𝔼_max(a)                = sum(μ⃗[jz] * max(ℑᴱ[jz](a), ℑᵂ[jz](a)) for jz in eachindex(z⃗))
        # ℑ(iz, a)                = β * (ψ * max(ℑᴱ[iz](a), ℑᵂ[iz](a)) + (1 - ψ) * 𝔼_max(a))
        @. mit_endo.𝔼𝐕[:,:,it]    = ψ * mit_endo.𝐕[:,:,it+1] + (1 - ψ) * $(μ⃗' * mit_endo.𝐕[:,:,it+1])

        # D2. Open the inner loops 
        Threads.@threads for ia in eachindex(a⃗)
            @inbounds for iz in eachindex(z⃗)

                # I. Find assets
                # ℑᶻ(a)   = ℑ(iz, a)
                # mit_endo.𝐚ʷ[iz, ia, it], mit_endo.𝐚ᵉ[iz, ia, it], mit_endo.𝐕ᵂ[iz, ia, it], mit_endo.𝐕ᴱ[iz, ia, it], 𝔼𝐕ᵂ[iz, ia, it], 𝔼𝐕ᴱ[iz, ia, it] = fnFindAssetsMIT(iz, ia, it, ℑᶻ, params, mit_endo)
                Yᵂ                                                  = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.wₜ[it] - mit_endo.τₜ[it]
                Yᴱ                                                  = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.Π[iz, ia,it] - mit_endo.τₜ[it]
                mit_endo.𝐚ʷ[iz, ia,it], mit_endo.𝐕ᵂ[iz, ia, it]     = fnPureGridSearchMIT(Yᵂ, params, a⃗, @views(mit_endo.𝔼𝐕[iz,:,it]))
                mit_endo.𝐚ᵉ[iz, ia,it], mit_endo.𝐕ᴱ[iz, ia, it]     = fnPureGridSearchMIT(Yᴱ, params, a⃗, @views(mit_endo.𝔼𝐕[iz,:,it]))


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
    B                           = zeros(Nᵃ, Nˡ, Nᵘ)
    gⁿᵉˣᵗ                       = zeros(Nᶻ, Nᵃ, Nˡ, Nᵘ)
    @. mit_endo.g[:,:,:,:, 1]   = @views(ss_endo.g)

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
                        Mass = mit_endo.g[iz, ia, il, iu, it]
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
                        B[ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 2]      += Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowU
                        B[ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 2]      += Mass * wᵃ[iz,ia, it] * (1 - wˡ[iz,ia,it]) * FlowU
                        B[iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 2]      += Mass * (1 - wᵃ[iz,ia,it]) * wˡ[iz,ia,it] * FlowU
                        B[iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 2]      += Mass * (1 - wᵃ[iz,ia,it]) * (1 - wˡ[iz,ia,it]) * FlowU
                        # To employment/entrepreneurship (iu=1)
                        B[ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 1]      += Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowO
                        B[ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 1]      += Mass * wᵃ[iz,ia,it] * (1 - wˡ[iz,ia,it]) * FlowO
                        B[iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 1]      += Mass * (1 - wᵃ[iz,ia,it]) * wˡ[iz,ia,it] * FlowO
                        B[iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 1]      += Mass * (1 - wᵃ[iz,ia,it]) * (1 - wˡ[iz,ia,it]) * FlowO

                        # V. The same productivity
                        # (c) Update (asset, labour) = (BB, BU, UB, UU)
                        # To unemployment 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 2]      += ψ * Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowU
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 2]      += ψ * Mass * wᵃ[iz,ia,it] * (1 - wˡ[iz,ia,it]) * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 2]      += ψ * Mass * (1 - wᵃ[iz,ia, it]) * wˡ[iz,ia, it] * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 2]      += ψ * Mass * (1 - wᵃ[iz,ia, it]) * (1 - wˡ[iz,ia,it]) * FlowU
                        # To others 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], ibˡ[iz,ia,it], 1]      += ψ * Mass * wᵃ[iz,ia,it] * wˡ[iz,ia,it] * FlowO
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia,it], iuˡ[iz,ia,it], 1]      += ψ * Mass * wᵃ[iz,ia,it] * (1 - wˡ[iz,ia,it]) * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], ibˡ[iz,ia,it], 1]      += ψ * Mass * (1 - wᵃ[iz,ia,it]) * wˡ[iz,ia,it] * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia,it], iuˡ[iz,ia,it], 1]      += ψ * Mass * (1 - wᵃ[iz,ia,it]) * (1 - wˡ[iz,ia,it]) * FlowO
                    end
                end
            end
        end

        # VI. Add the mass to the updated productivity states
        @inbounds for iu in 1:2, il in 1:Nˡ ,ia in 1:Nᵃ, iz in 1:Nᶻ
            gⁿᵉˣᵗ[iz, ia, il, iu]   += (1 - ψ) * μ⃗[iz] * B[ia, il, iu]
        end

        # VII. Update 
        @. mit_endo.g[:,:,:,:,it+1] = gⁿᵉˣᵗ
        fill!(gⁿᵉˣᵗ, 0.0)
        fill!(B, 0.0)
    end

    # Final update to the important rates 
    fnJobDestructionMIT!(params, mit_endo, Tᴹᴵᵀ)
    fnUpdateLabourMarketMIT!(params, mit_endo, Tᴹᴵᵀ)
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

# -----------------------------------------------
# C. Solving
# -----------------------------------------------

# 1. Convex updating for MIT shock
function fnConvexUpdatingMIT(f, bounds::Tuple, method = nothing; 
                            loading = 0.05, xatol = 1e-5, max_iter=2000, 
                            init = nothing, step_tol = 1e-3, kwargs...)
    
    # A. Initialization
    xˣ = (isnothing(init) || all(init .== 0.0)) ? 
         fill((bounds[1] + bounds[2]) / 2.0, length(bounds[1])) : copy(init)
    
    loadings = fill(loading, length(xˣ))
    prev_res = zeros(length(xˣ))
    
    iter    = 0
    ε       = 1.0
    x_old   = fill(0.0, length(xˣ))

    # B. Start the loop 
    while ε > xatol && iter < max_iter
        residual        = f(xˣ)
        
        # C. Limit 
        if iter > 0
            @. loadings = ifelse(sign(residual) != sign(prev_res), loadings * 0.7, loadings)
        end

        # D. Safety clause
        if maximum(loadings) < loading/10 
            return xˣ 
        end
        if iter>3 && maximum(abs, xˣ .- x_old) < step_tol
            @warn "Solver stagnated at step $iter. Exiting."
            return xˣ
        end

        # E. Update 
        x_old       .= copy(xˣ)
        @. xˣ       = clamp(xˣ * (1 + loadings * residual), bounds[1], bounds[2])
        prev_res    .= residual
        ε           = maximum(abs, residual)
        iter        += 1
    end

    if iter == max_iter
        @warn "MIT Updating reached max_iter. ε = $ε"
    end

    return xˣ
end

# 2. Labour market residual (MIT)
function fnLabourResidualMIT(w⃗,r⃗, τ⃗, params, mit_endo, ss_endo,  A⃗, λ⃗)
    
    # A. Unpacking parameters
    @unpack Tᴹᴵᵀ = params

    # B. Update the state
    mit_endo.wₜ = w⃗
    mit_endo.rₜ = r⃗
    mit_endo.τₜ = τ⃗

    # C. Run the VFI and aggregate routine 
    fnBackwardInductionMIT!(params, mit_endo, ss_endo, A⃗, λ⃗)
    fnAggregateStatesMIT!(params, mit_endo, ss_endo)
    
    # D. Return the labour market error
    ε⃗ᴸ          = @. (mit_endo.Lᵈ / max(mit_endo.Lˢ, 1e-4)) - 1.0
    εᴸ          = maximum(ε⃗ᴸ) 
    ε̃ᴸ          = minimum(ε⃗ᴸ) 
    𝔼w          = sum(mit_endo.wₜ) / Tᴹᴵᵀ
    println("→→→ Wage search        | w₄ = $(round(mit_endo.wₜ[4], digits=4)), 𝔼w = $(round(𝔼w, digits=4)), max(εᴸ) = $(round(εᴸ, digits=4)), min(εᴸ) = $(round(ε̃ᴸ, digits=4))")
    return ε⃗ᴸ
end 

# 2. Budget residual (MIT)
# A. A struct to hold the variables needed for the labor residual
struct LabourResidualObjectiveMIT{P,T,E1,E2,G}
    params::P
    r⃗::T
    τ⃗::T
    mit_endo::E1
    ss_endo::E2
    A⃗::G 
    λ⃗::G
end

# B. Make it callable
function (obj::LabourResidualObjectiveMIT)(w⃗)
    return fnLabourResidualMIT(w⃗, obj.r⃗, obj.τ⃗, obj.params, obj.mit_endo,obj.ss_endo,obj.A⃗,obj.λ⃗)
end

# C. Start the function 
function fnBudgetResidualMIT!(τ⃗,r⃗, params, mit_endo,ss_endo,A⃗, λ⃗)

    # D. Unpacking business 
    @unpack δᴸ, κᴸ,Tᴹᴵᵀ = params

    # E. Prepare lower and upper bounds 
    w̲       = fill(0.1,Tᴹᴵᵀ)
    w̅       = fill(2.5,Tᴹᴵᵀ)

    # F. Find the wage that clears labor for THIS tax and interest rate
    Oᴸ              = LabourResidualObjectiveMIT(params, r⃗, τ⃗, mit_endo,ss_endo,A⃗, λ⃗)
    wˣ              = fnConvexUpdatingMIT(Oᴸ, (w̲, w̅),loading = κᴸ,xatol = 10*δᴸ,init = mit_endo.wₜ,max_iter=20)
    mit_endo.wₜ     .= wˣ 

    # G. Return the budget error
    ε⃗ᵗ          = @. (mit_endo.U * mit_endo.wₜ / max(mit_endo.τₜ, 1e-4)) - 1.0
    εᵗ          = maximum(ε⃗ᵗ)
    ε̃ᵗ          = minimum(ε⃗ᵗ)
    𝔼τ          = sum(τ⃗) / Tᴹᴵᵀ
    println("→→ Tax loop            | τ₄ = $(round(τ⃗[4], digits=4)), 𝔼τ = $(round(𝔼τ, digits=4)), max(εᵗ) = $(round(εᵗ, digits=4)),  min(εᵗ) = $(round(ε̃ᵗ, digits=4))")
    return ε⃗ᵗ
end

# 3. Government budget 
# A. Struct for the tax residual
struct BudgetResidualObjectiveMIT{P,T,E1,E2,G}
    params::P
    r⃗::T
    mit_endo::E1
    ss_endo::E2
    A⃗::G 
    λ⃗::G
end

# B. Make it callable
function (obj::BudgetResidualObjectiveMIT)(τ⃗)
    return fnBudgetResidualMIT!(τ⃗, obj.r⃗,obj.params, obj.mit_endo,obj.ss_endo,obj.A⃗,obj.λ⃗)
end

function fnCapitalResidualMIT!(r⃗, params, mit_endo,ss_endo,A⃗, λ⃗, error_history)
        
        # C. Unpacking business 
        @unpack δᵗ,κᵗ,Tᴹᴵᵀ,δᴸ,κᴸ = params

        # D. Prepare lower and upper bounds 
        τ̲       = fill(0.01,Tᴹᴵᵀ)
        τ̅       = fill(0.3,Tᴹᴵᵀ)
        
        # E. Solve [old, super precise but at risk of looping]
        Oᵗ      = BudgetResidualObjectiveMIT(params, r⃗, mit_endo,ss_endo,A⃗, λ⃗)
        τˣ      = fnConvexUpdatingMIT(Oᵗ, (τ̲, τ̅),loading = 0.5*κᵗ,  xatol = 10*δᵗ,init = mit_endo.τₜ)
        @. mit_endo.τₜ      = τˣ

        # E. Solve: iterate w-search → τ update a few times
        # for _ in 1:3
        #     Oᴸ              = LabourResidualObjectiveMIT(params, r⃗, mit_endo.τₜ, mit_endo, ss_endo, A⃗, λ⃗)
        #     wˣ              = fnConvexUpdatingMIT(Oᴸ, (fill(0.1, Tᴹᴵᵀ), fill(2.5, Tᴹᴵᵀ)), loading=κᴸ, xatol=δᴸ, init=mit_endo.wₜ)
        #     mit_endo.wₜ     .= wˣ
        #     @. mit_endo.τₜ  = mit_endo.wₜ * mit_endo.U
        # end

        # F. Return error 
        ε⃗ᴷ  = @. (mit_endo.Kᵈ / max(mit_endo.Kˢ, 1e-4)) - 1.0
        εᴷ  = maximum(abs,ε⃗ᴷ)
        push!(error_history, εᴷ)
        fnPlotConvergenceMIT(error_history, mit_endo, ss_endo,params)
        println("⋆ Capital loop     | r₄ = $(round(r⃗[4], digits=4)), max(εᴷ) = $(round(εᴷ, digits=4))")
        return ε⃗ᴷ
end

# 4. Steady state 

# A. Struct for the capital residual
struct CapitalResidualObjectiveMIT{P,E1,E2,G,H}
    params::P
    mit_endo::E1
    ss_endo::E2
    A⃗::G 
    λ⃗::G
    error_history::H
end

# B. Make it callable
function (obj::CapitalResidualObjectiveMIT)(r⃗)
    return fnCapitalResidualMIT!(r⃗, obj.params, obj.mit_endo,obj.ss_endo,obj.A⃗,obj.λ⃗,obj.error_history)
end

# C. Start the function 
function fnTransitionMIT!(params, mit_endo, ss_endo, A⃗, λ⃗)
    
    # D. Unpacking business 
    @unpack δʳ,κʳᴹᴵᵀ,Tᴹᴵᵀ,δ,r̲,r̅ = params

    # E. Bounds and warm start 
    r̲⃗               = fill(-δ,Tᴹᴵᵀ)
    r̅⃗               = fill(0.20,Tᴹᴵᵀ)
    # fill!(mit_endo.rₜ, ss_endo.rₜ)
    # fill!(mit_endo.wₜ, ss_endo.wₜ)
    Multiplier    = min.(A⃗ ./ A⃗[1], λ⃗ ./ λ⃗[1])
    fill!(mit_endo.τₜ, ss_endo.τₜ)
    mit_endo.rₜ  .= ss_endo.rₜ .* Multiplier
    mit_endo.wₜ  .= ss_endo.wₜ .* Multiplier

    # F. The final solve
    error_history       = Float64[]
    Oᶜ                  = CapitalResidualObjectiveMIT(params, mit_endo, ss_endo, A⃗, λ⃗,error_history)
    # logit(x, a, b)      = log((x - a) / (b - x))
    # invlogit(y, a, b)   = a + (b - a) / (1 + exp(-y))
    # y₀                  = logit.(mit_endo.rₜ, r̲⃗, r̅⃗)
    # result              = nlsolve((F, y) -> F .= Oᶜ(invlogit.(y, r̲⃗, r̅⃗)), y₀, method=:anderson, m=5, ftol=10*δʳ)
    # mit_endo.rₜ         .= invlogit.(result.zero, r̲⃗, r̅⃗)

    rˣ              = fnConvexUpdatingMIT(Oᶜ, (r̲, r̅),loading=κʳᴹᴵᵀ, xatol = δʳ, init = mit_endo.rₜ,step_tol=1e-4)
    @. mit_endo.rₜ  = rˣ
    # G. Save convergence plot
    plt_conv = fnPlotConvergenceMIT(error_history, mit_endo,ss_endo, params)
    savefig(plt_conv, "plots/MIT_convergence.pdf")
end

# 5. Live plotting function for MIT transition
function fnPlotConvergenceMIT(error_history, mit_endo, ss_endo, params)
    
    @unpack Tᴹᴵᵀ = params
    t_grid = 1:Tᴹᴵᵀ

    # Row 1
    p1 = plot(error_history, title="Capital error history", xlabel="t", ylabel="Max absolute error", yscale=:log10, lw=2, color=:black, legend=false, grid=true)
    
    p2 = plot(t_grid, mit_endo.Kᵈ, title="Capital", xlabel="t", lw=2, color=:black, label="Demand", grid=true)
    plot!(p2, t_grid, mit_endo.Kˢ, lw=2, color=:red, ls=:dash, label="Supply")
    hline!(p2, [ss_endo.Kˢ], color=:gray, ls=:dash, label="Steady state")
    
    p3 = plot(t_grid, mit_endo.rₜ, title="Interest rate", xlabel="t", lw=2, color=:blue, legend=false, grid=true)
    hline!(p3, [ss_endo.rₜ], color=:gray, ls=:dash)

    # Row 2
    p4 = plot(t_grid, mit_endo.Lᵈ, title="Labour", xlabel="t", lw=2, color=:black, label="Demand", grid=true)
    plot!(p4, t_grid, mit_endo.Lˢ, lw=2, color=:red, ls=:dash, label="Supply")
    hline!(p4, [ss_endo.Lˢ], color=:gray, ls=:dash, label="SS")

    p5 = plot(t_grid, mit_endo.U, title="Unemployment", xlabel="t", lw=2, color=:orange, legend=false, grid=true)
    hline!(p5, [ss_endo.U], color=:gray, ls=:dash)

    p6 = plot(t_grid, mit_endo.wₜ, title="Wage", xlabel="t", lw=2, color=:green, legend=false, grid=true)
    hline!(p6, [ss_endo.wₜ], color=:gray, ls=:dash)

    # Row 3
    p7 = plot(t_grid, mit_endo.τₜ, title="Tax", xlabel="t", lw=2, color=:red, legend=false, grid=true)
    hline!(p7, [ss_endo.τₜ], color=:gray, ls=:dash)

    p8 = plot(t_grid, mit_endo.E, title="Entrepreneurs", xlabel="t", lw=2, color=:purple, legend=false, grid=true)
    hline!(p8, [ss_endo.E], color=:gray, ls=:dash)

    p9 = plot(t_grid, mit_endo.JD, title="Job Destruction", xlabel="t", lw=2, color=:black, label="Total", grid=true)
    plot!(p9, t_grid, mit_endo.S, lw=2, color=:blue, ls=:dash, label="Switchers (S)")
    plot!(p9, t_grid, mit_endo.D, lw=2, color=:red, ls=:dot, label="Downsizing (D)")
    hline!(p9, [ss_endo.JD], color=:gray, ls=:dash, label="SS Total")

    plt = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(3, 3), size=(1200, 1000), margin=5Plots.mm)
    display(plt)
    return plt
end