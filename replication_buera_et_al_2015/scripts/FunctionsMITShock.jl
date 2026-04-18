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
 function fnPureGridSearchMIT(Budget, params, a⃗, 𝔼V_array,LowerBound_i)
    
    # A. Unpacking business
    @unpack c̲, β, Nᵃ = params
    
    # B. Initial values 
    best_val    = -Inf
    best_a      = a⃗[LowerBound_i]
    best_i      = LowerBound_i
    
    # C. Start the loop 
    @inbounds for i in LowerBound_i:Nᵃ
        
        # D. Savings and consumption 
        a_prime         = a⃗[i]
        c               = Budget - a_prime
        
        # D. Stop searching if consumption hits the minimum bound
        if c <= c̲
            v_bound             = fnUtility(c̲, params) + β * 𝔼V_array[i]
            if v_bound > best_val
                best_val        = v_bound
                best_a          = a_prime
                best_i          = LowerBound_i
            end
            break 
        end
        
        # E. Pure array lookup 
        v                       = fnUtility(c, params) + β * 𝔼V_array[i]
        if v > best_val
            best_val            = v
            best_a              = a_prime
            best_i              = i
        end
    end
    
    # F. Returning business 
    return best_a, best_val, best_i
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
        μᵥ                      = μ⃗' * mit_endo.𝐕[:,:,it+1]   
        @. mit_endo.𝔼𝐕[:,:,it]  = ψ * mit_endo.𝐕[:,:,it+1] + (1 - ψ) * μᵥ

        # D2. Open the inner loops 
        Threads.@threads for iz in eachindex(z⃗)
            iʷ  = 1
            iᵉ  = 1

            @inbounds for ia in eachindex(a⃗)

                # I. Find assets
                # ℑᶻ(a)   = ℑ(iz, a)
                # mit_endo.𝐚ʷ[iz, ia, it], mit_endo.𝐚ᵉ[iz, ia, it], mit_endo.𝐕ᵂ[iz, ia, it], mit_endo.𝐕ᴱ[iz, ia, it], 𝔼𝐕ᵂ[iz, ia, it], 𝔼𝐕ᴱ[iz, ia, it] = fnFindAssetsMIT(iz, ia, it, ℑᶻ, params, mit_endo)
                Yᵂ                                                  = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.wₜ[it] - mit_endo.τₜ[it]
                Yᴱ                                                  = a⃗[ia] * (1 + mit_endo.rₜ[it]) + mit_endo.Π[iz, ia,it] - mit_endo.τₜ[it]
                mit_endo.𝐚ʷ[iz, ia,it], mit_endo.𝐕ᵂ[iz, ia, it],iʷ₊ = fnPureGridSearchMIT(Yᵂ, params, a⃗, @views(mit_endo.𝔼𝐕[iz,:,it]),iʷ)
                mit_endo.𝐚ᵉ[iz, ia,it], mit_endo.𝐕ᴱ[iz, ia, it],iᵉ₊ = fnPureGridSearchMIT(Yᴱ, params, a⃗, @views(mit_endo.𝔼𝐕[iz,:,it]),iᵉ)
                iʷ  = iʷ₊
                iᵉ  = iᵉ₊

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
    WU  = 0.0 + (il== 1) * JDR * (1 - JFR) * (iu == 1) * (mit_endo.𝐨[iz, ia, it] == false)
    WE  = 0.0 + (il== 1) * (mit_endo.𝐨[iz, ia, it] == true) * (iu == 1)
    EU  = 0.0 + (il > 1) * (mit_endo.𝐨[iz, ia, it] == false) * (iu == 1) * (1 - JFR)
    EW  = 0.0 + (il > 1) * (mit_endo.𝐨[iz, ia, it] == false) * (iu == 1) * JFR
    UE  = 0.0 + (mit_endo.𝐨[iz, ia, it] == true) * (iu == 2)
    UW  = 0.0 + (mit_endo.𝐨[iz, ia, it] == false) * (iu == 2) * JFR
    EE  = 0.0 + (il > 1) * (mit_endo.𝐨[iz, ia, it] == true) * (iu == 1)
    WW  = 0.0 + (il== 1) * (mit_endo.𝐨[iz, ia, it] == false) * (iu == 1) * ((1 - JDR) + (JDR * JFR))
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
    
    # D. Start the largest loop  ibₐ, iuₐ, etc
    @inbounds for it in 1:Tᴹᴵᵀ-1

        # I. Update the endogenous functions for time t
        fnJobDestructionMIT!(params, mit_endo, it)
        fnUpdateLabourMarketMIT!(params, mit_endo, it)

        for iu in 1:Nᵘ
            for il in eachindex(l⃗)
                for ia in eachindex(a⃗)
                    for iz in eachindex(z⃗)

                        # II. Common terms
                        Mass            = mit_endo.g[iz, ia, il, iu, it]
                        Mass            == 0.0 && continue

                        # III. Hoist indices and weights once
                        ibₐ, iuₐ, wₐ    = ibᵃ[iz,ia,it], iuᵃ[iz,ia,it], wᵃ[iz,ia,it]
                        ibₗ, iuₗ, wₗ    = ibˡ[iz,ia,it], iuˡ[iz,ia,it], wˡ[iz,ia,it]

                        # IV. Flows
                        f¹², f²¹, f²², f¹¹  = fnComputeFlowsMIT(iz, ia, il, iu, it, params, mit_endo)
                        FlowU, FlowO        = f¹² + f²², f²¹ + f¹¹

                        # V. Precomputed corner weights (4 mults instead of 16)
                        wbb     = wₐ * wₗ
                        wbu     = wₐ * (1 - wₗ)
                        wub     = (1 - wₐ) * wₗ
                        wuu     = (1 - wₐ) * (1 - wₗ)
                        MU, MO = Mass * FlowU, Mass * FlowO

                        # VI. Write to B (switching productivity block)
                        B[ibₐ, ibₗ, 2]  += MU * wbb;  B[ibₐ, iuₗ, 2] += MU * wbu
                        B[iuₐ, ibₗ, 2]  += MU * wub;  B[iuₐ, iuₗ, 2] += MU * wuu
                        B[ibₐ, ibₗ, 1]  += MO * wbb;  B[ibₐ, iuₗ, 1] += MO * wbu
                        B[iuₐ, ibₗ, 1]  += MO * wub;  B[iuₐ, iuₗ, 1] += MO * wuu

                        # VII. gⁿᵉˣᵗ gets ψ × same contributions at own iz slice
                        ψMU, ψMO = ψ * MU, ψ * MO
                        gⁿᵉˣᵗ[iz, ibₐ, ibₗ, 2]  += ψMU * wbb;  gⁿᵉˣᵗ[iz, ibₐ, iuₗ, 2] += ψMU * wbu
                        gⁿᵉˣᵗ[iz, iuₐ, ibₗ, 2]  += ψMU * wub;  gⁿᵉˣᵗ[iz, iuₐ, iuₗ, 2] += ψMU * wuu
                        gⁿᵉˣᵗ[iz, ibₐ, ibₗ, 1]  += ψMO * wbb;  gⁿᵉˣᵗ[iz, ibₐ, iuₗ, 1] += ψMO * wbu
                        gⁿᵉˣᵗ[iz, iuₐ, ibₗ, 1]  += ψMO * wub;  gⁿᵉˣᵗ[iz, iuₐ, iuₗ, 1] += ψMO * wuu
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

# 1. Implied labour supply error 
function fnImpliedLabourDemandError(w, params, mit_endo, A⃗, λ⃗, it)
    # A. Unpacking
    @unpack α, δ, θ, z⃗, a⃗ = params
    
    # B. Individual values 
    A1 = @. α * A⃗[it] * z⃗ / (mit_endo.rₜ[it] + δ)
    A2 = θ / α * (mit_endo.rₜ[it] + δ) / w   
    A3 = A2^θ
    kˣ   = @. (A1 * A3)^(1 / (1 - α - θ))
    kⁱᵐᵖ = @. min(kˣ, λ⃗[it] * a⃗')
    lⁱᵐᵖ = @. (θ * A⃗[it] * z⃗ * kⁱᵐᵖ^α / w)^(1 / (1 - θ)) # 

    # C. Aggregate demand
    ωᵃ   = dropdims(sum(@views(mit_endo.g[:, :, :, :, it]), dims=(3, 4)), dims=(3, 4))
    Lⁱᵐᵖ = sum(ωᵃ .* lⁱᵐᵖ .* @views(mit_endo.𝐨[:, :, it]))

    # D. Return error (Supply - Demand for this specific period)
    return mit_endo.Lˢ[it] - Lⁱᵐᵖ  
end

# 2. Implied capital demand error 
function fnImpliedCapitalDemandError(r, params, mit_endo, A⃗, λ⃗, it)
    # A. Unpacking
    @unpack α, δ, θ, z⃗, a⃗ = params
    
    # B. Individual values 
    A1 = @. α * A⃗[it] * z⃗ / (r + δ)          
    A2 = θ / α * (r + δ) / mit_endo.wₜ[it]   
    A3 = A2^θ
    kˣ   = @. (A1 * A3)^(1 / (1 - α - θ))
    kⁱᵐᵖ = @. min(kˣ, λ⃗[it] * a⃗')            

    # C. Aggregate demand
    ωᵃ   = dropdims(sum(@views(mit_endo.g[:, :, :, :, it]), dims=(3, 4)), dims=(3, 4))
    Kⁱᵐᵖ = sum(ωᵃ .* kⁱᵐᵖ .* @views(mit_endo.𝐨[:, :, it]))

    # D. Return error 
    return mit_endo.Kˢ[it] - Kⁱᵐᵖ
end

# 3. MIT solver 
function fnSolveMIT!(params, mit_endo, ss_endo, A⃗, λ⃗)

    # A. Unpacking business 
    @unpack Tᴹᴵᵀ, ηˡ, ηᵗ, ηᶜ, δᴸ, δᵗ, δʳ, δ, β = params 

    # B. Initial guesses 
    error_history   = Float64[]
    for it in 1:Tᴹᴵᵀ
        mit_endo.rₜ[it] = ss_endo.rₜ * (λ⃗[it] / λ⃗[1])^0.5
    end 
    mit_endo.wₜ     .= ss_endo.wₜ 
    mit_endo.τₜ     .= mit_endo.wₜ .* ss_endo.U

    # C. Interest and wage limits (for bisections)
    w̲       = 0.4
    w̅       = 2.5
    r̲       = -δ +1e-6
    r̅       = 0.15

    # D. Loop settings 
    εʳ      = 1.0
    εᵗ      = 1.0
    εˡ      = 1.0  
    w̃       = zeros(Tᴹᴵᵀ)
    τ̃       = zeros(Tᴹᴵᵀ)
    r̃       = zeros(Tᴹᴵᵀ)
    # Temporary values 
    δᴸ      = 5*1e-4
    δᵗ      = 1e-3
    δʳ      = 5*1e-3

    # E1. Starting the loops 
    while εʳ > δʳ
        εᵗ  = 1.0 
        η⃗ˡ  = fill(ηˡ, Tᴹᴵᵀ)
        𝓃ᵗ  = 0
        𝓃̄ᵗ  = 15

        # E2. Start the tax loop 
        while εᵗ > δᵗ && 𝓃ᵗ <= 𝓃̄ᵗ
            εˡ  = 1.0
            𝓃ʷ  = 0
            𝓃̄ʷ = εʳ > 0.05 ? 50 : 20

            # E3. Start the wage loop 
            while εˡ > δᴸ && 𝓃ʷ < 𝓃̄ʷ

                # I. Solve the model 
                tᵇ = @elapsed fnBackwardInductionMIT!(params, mit_endo, ss_endo, A⃗, λ⃗)
                tᶠ = @elapsed fnAggregateStatesMIT!(params, mit_endo, ss_endo)

                # II. Implied wage 
                for it in 1:Tᴹᴵᵀ
                    flo = fnImpliedLabourDemandError(w̲, params, mit_endo, A⃗, λ⃗, it)
                    fhi = fnImpliedLabourDemandError(w̅, params, mit_endo, A⃗, λ⃗, it)
                    if flo * fhi > 0
                        @warn "No bracket at t=$it: f(w̲)=$flo, f(w̅)=$fhi"
                    end
                    w̃[it] = find_zero(w -> fnImpliedLabourDemandError(w, params, mit_endo, A⃗, λ⃗, it), (w̲, w̅), Bisection())
                end 

                # III. Error term 
                ε⃗ᴸ      = @. (mit_endo.Lᵈ / max(mit_endo.Lˢ, 1e-4)) - 1.0
                εˡ      = maximum(abs, ε⃗ᴸ) 
                
                # IV. Update
                @. η⃗ˡ = ifelse(abs(ε⃗ᴸ) < 0.0025, 0.25, ifelse(abs(ε⃗ᴸ) < 0.025, 0.15, ifelse(abs(ε⃗ᴸ) < 0.075, 0.5, ηˡ)))
                @. mit_endo.wₜ = η⃗ˡ * w̃ + (1 - η⃗ˡ) * mit_endo.wₜ

                # V. Print a fun message for future generations 
                println("→→→ Wage search        | w₄ = $(round(mit_endo.wₜ[4], digits=4)), max(εᴸ) = $(round(εˡ, digits=7)), tᵇ = $(round(tᵇ,digits=3)), tᶠ = $(round(tᶠ,digits=3))")
            
                # VI. Iteration limits 
                𝓃ʷ += 1
                if 𝓃ʷ == 𝓃̄ʷ
                    println("⚠ Wage cap hit: εᴸ = $εˡ")
                end
            end 

            # VI. Calculate the implied tax 
            @. τ̃            = mit_endo.U * mit_endo.wₜ

            # VII. Error term 
            ε⃗ᵗ              = @. (mit_endo.U * mit_endo.wₜ / max(mit_endo.τₜ, 1e-4)) - 1.0
            εᵗ              = maximum(abs, ε⃗ᵗ) 

            # VIII. Update 
            @. mit_endo.τₜ  = ηᵗ * τ̃ + (1 - ηᵗ) * mit_endo.τₜ

            # IX. Nice message 
            println("→→ Tax loop            | τ₄ = $(round(mit_endo.τₜ[4], digits=4)), max(εᵗ) = $(round(εᵗ, digits=6))")
        end 

        # X. Implied interest rate 
        for it in 1:Tᴹᴵᵀ
            try
                r̃[it] = find_zero(r -> fnImpliedCapitalDemandError(r, params, mit_endo, A⃗, λ⃗, it), (r̲, r̅), Bisection())
            catch
                flo = fnImpliedCapitalDemandError(r̲, params, mit_endo, A⃗, λ⃗, it)
                fhi = fnImpliedCapitalDemandError(r̅, params, mit_endo, A⃗, λ⃗, it)
                @show flo, fhi, r̲, r̅
                r̃[it] = r̲  # keep current guess
                @warn "No bracket at t=$it, holding r"
            end
        end

        # XI. Error term and plot update 
        ε⃗ᴷ  = @. (mit_endo.Kᵈ / max(mit_endo.Kˢ, 1e-4)) - 1.0
        εʳ  = maximum(abs,ε⃗ᴷ)
        push!(error_history, εʳ)
        fnPlotConvergenceMIT(error_history, mit_endo, ss_endo,params,r̃)

        # XII. Update 
        @. mit_endo.rₜ  = ηᶜ * r̃ + (1 - ηᶜ) * mit_endo.rₜ

        # XIII. Lovely message
        println("⋆ Capital loop     | r₄ = $(round(mit_endo.rₜ[4], digits=4)), max(εᴷ) = $(round(εʳ, digits=6))") 
    end
end 

# 4. Live plotting function for MIT transition
function fnPlotConvergenceMIT(error_history, mit_endo, ss_endo, params, r_implied)
    
    @unpack Tᴹᴵᵀ = params
    t_grid = 1:Tᴹᴵᵀ

    # Row 1
    p1 = plot(error_history, title="Capital error history", xlabel="t", ylabel="Max absolute error", yscale=:log10, lw=2, color=:black, legend=false, grid=true)
    
    p2 = plot(t_grid, mit_endo.Kᵈ, title="Capital", xlabel="t", lw=2, color=:black, label="Demand", grid=true)
    plot!(p2, t_grid, mit_endo.Kˢ, lw=2, color=:red, ls=:dash, label="Supply")
    hline!(p2, [ss_endo.Kˢ], color=:gray, ls=:dash, label="Steady state")
    
    p3 = plot(t_grid, mit_endo.rₜ, title="Interest rate", xlabel="t", lw=2, color=:blue, legend=false, grid=true, label = "Guess")
    plot!(p3, t_grid, r_implied, lw=2, color=:blue, ls=:dash, label="Implied")
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
