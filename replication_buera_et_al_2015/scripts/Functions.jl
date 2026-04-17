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

# C. Solving  
# 1A. Logistic function for bounding 
# 1B. Inverse logistic for bounding
# 2. Combined steady state sesiduals
# 3. Solve for the steady state 

# -----------------------------------------------
# A. VFI components 
# -----------------------------------------------

# 1. Utility function 
function fnUtility(c, params)
    @unpack σ, c̲ = params
    return c <= c̲ ? (c̲^(1 - σ)) / (1 - σ) + (c - c̲) * (c̲^(-σ)) : (c^(1 - σ)) / (1 - σ)
end

# 2. Static policy functions  
function fnStaticPolicies!(params, endo)
    @unpack α, δ, θ, A, z⃗, λ, a⃗ = params

    A1 = @. α * A * z⃗ / (endo.rₜ + δ)
    A2 = θ / α * (endo.rₜ + δ) / endo.wₜ
    A3 = A2^(θ)

    # Cleaned up manual dots
    kˣ          = @. (A1 * A3) ^ (1 / (1 - α - θ))
    @. endo.𝐤   = min(kˣ, λ * a⃗')
    @. endo.𝕀ᶜ  = kˣ >= λ * a⃗'

    @. endo.𝐥   = (θ * A * z⃗ * endo.𝐤 ^ α / endo.wₜ) ^ (1 / (1 - θ))
    @. endo.Π   = A * z⃗ * endo.𝐤 ^ α * endo.𝐥 ^ θ - endo.wₜ * endo.𝐥 - (endo.rₜ + δ) * endo.𝐤
end

# 3. Initial VFI guesses 
function fnInitialVFIGuess!(params, endo)
    @unpack β, a⃗ = params

    if endo.𝐕[1, 1] == 0.0
        endo.𝐕ᵂ     .=  (1 - β)^(-1) .* fnUtility.(endo.rₜ .* a⃗' .- endo.τₜ .+ endo.wₜ, Ref(params))
        endo.𝐕ᴱ     .=  (1 - β)^(-1) .* fnUtility.(endo.rₜ .* a⃗' .- endo.τₜ .+ endo.Π, Ref(params))

        @. endo.𝐕   = max(endo.𝐕ᵂ, endo.𝐕ᴱ)
        @. endo.𝐨   = endo.𝐕ᵂ <= endo.𝐕ᴱ
    end 
end

# 4A. Solve for assets 
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
function fnFindAssets(iz, ia, exp_spline, params, endo)

    # D. Unpacking parameters
    @unpack c̲, a⃗, β = params

    # E. Helper elements 
    Cash_w = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ
    Cash_e = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz, ia] - endo.τₜ
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

# 4B. Optimised grid search 
 function fnPureGridSearch(Budget, params, a⃗, 𝔼V_array)
    
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

# 5. Value function iteration 
function fnVFI!(params, endo)

    # A. Unpacking business
    @unpack ψ, μ⃗, δᵛᶠⁱ, 𝒾̄ᵛᶠⁱ, z⃗, c̲, a⃗, β, λᵛᶠⁱ, Nᶻ, Nᵃ, interp = params
    
    # B. Static policies
    fnStaticPolicies!(params, endo)

    # C. Initial guess 
    fnInitialVFIGuess!(params, endo)

    # D. Prepare the VFI loop 
    εᵛᶠⁱ    = 1.0
    𝓃ᵛᶠⁱ    = 1
    𝔼𝐕ᵂ     = zeros(Nᶻ, Nᵃ) 
    𝔼𝐕ᴱ     = zeros(Nᶻ, Nᵃ)
    𝐕ᵖʳᵉᵛ   = copy(endo.𝐕)

    while (εᵛᶠⁱ > δᵛᶠⁱ && 𝓃ᵛᶠⁱ < 𝒾̄ᵛᶠⁱ)
        𝐕ᵖʳᵉᵛ       = copy(endo.𝐕)
        if interp == true
            # D1. Make safe, static copies of the value functions for the splines to read
            𝐕ᴱ_read     = copy(endo.𝐕ᴱ)
            𝐕ᵂ_read     = copy(endo.𝐕ᵂ)

            # D2. Interpolate the expected value properly using the safe copies.
            ℑᴱ          = [linear_interpolation(a⃗, 𝐕ᴱ_read[iz, :], extrapolation_bc=Interpolations.Flat()) for iz in eachindex(z⃗)]
            ℑᵂ          = [linear_interpolation(a⃗, 𝐕ᵂ_read[iz, :], extrapolation_bc=Interpolations.Flat()) for iz in eachindex(z⃗)]

            # D3. Define the z-independent expectation closure outside the threaded loop.
            𝔼_max(a)    = sum(μ⃗[jz] * max(ℑᴱ[jz](a), ℑᵂ[jz](a)) for jz in eachindex(z⃗))

            # D4. The final objective closure passed to the optimizer
            ℑ(iz, a)    = β * (ψ * max(ℑᴱ[iz](a), ℑᵂ[iz](a)) + (1 - ψ) * 𝔼_max(a))

            # D5. Update the expected VF and open the loop for productivity 
        else
            μᵥ          = μ⃗' * endo.𝐕
            @. endo.𝔼𝐕  = ψ * endo.𝐕 + (1 - ψ) * μᵥ
        end

        Threads.@threads for ia in eachindex(a⃗)
            @inbounds for iz in eachindex(z⃗)
                # B1. Solve for each a and z 
                # I. Find assets
                if interp == true
                    let iz⁰ = iz
                        ℑᶻ  = a -> ℑ(iz⁰, a)
                        endo.𝐚ʷ[iz, ia], endo.𝐚ᵉ[iz, ia], endo.𝐕ᵂ[iz, ia], endo.𝐕ᴱ[iz, ia], 𝔼𝐕ᵂ[iz, ia], 𝔼𝐕ᴱ[iz, ia] = fnFindAssets(iz, ia, ℑᶻ, params, endo)
                    end
                else
                    # Budget for workers and entrepreneurs 
                    Yᵂ = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ
                    Yᴱ = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz, ia] - endo.τₜ
                    endo.𝐚ʷ[iz, ia], endo.𝐕ᵂ[iz, ia] = fnPureGridSearch(Yᵂ, params, a⃗, @views(endo.𝔼𝐕[iz, :]))
                    endo.𝐚ᵉ[iz, ia], endo.𝐕ᴱ[iz, ia] = fnPureGridSearch(Yᴱ, params, a⃗, @views(endo.𝔼𝐕[iz, :]))
                end

                # II. VFs for different occupations
                endo.𝐜ʷ[iz, ia] = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ - endo.𝐚ʷ[iz, ia]
                endo.𝐜ᵉ[iz, ia] = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz, ia] - endo.τₜ - endo.𝐚ᵉ[iz, ia]

                # III. Occupational decision and updated VF
                endo.𝐨[iz, ia] = (endo.𝐕ᴱ[iz, ia] >= endo.𝐕ᵂ[iz, ia])
                endo.𝐚[iz, ia] = endo.𝐨[iz, ia] ? endo.𝐚ᵉ[iz, ia] : endo.𝐚ʷ[iz, ia]
                endo.𝐜[iz, ia] = max(endo.𝐨[iz, ia] ? endo.𝐜ᵉ[iz, ia] : endo.𝐜ʷ[iz, ia], c̲)
                endo.𝐕[iz, ia] = max(endo.𝐕ᴱ[iz, ia], endo.𝐕ᵂ[iz, ia])
            end
        end

        # D6. Error and update 
        εᵛᶠⁱ = maximum(abs.(𝐕ᵖʳᵉᵛ .- endo.𝐕))
        𝓃ᵛᶠⁱ += 1
        if 𝓃ᵛᶠⁱ % 100 == 0
            println("→→→→ VFI Iteration: $𝓃ᵛᶠⁱ, εᵛᶠⁱ = $(round(εᵛᶠⁱ, digits=8))")
        end
        if 𝓃ᵛᶠⁱ == 𝒾̄ᵛᶠⁱ
            println("Warning: VFI reached the maximum number of iterations")
        end
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
    endo.S = 0.0
    endo.D = 0.0

    # C. Loop 
    for iz in eachindex(z⃗)
        for ia in eachindex(a⃗)
            for il in eachindex(l⃗)

                # C1. Compute key elements 
                Mass            = endo.g[iz, ia, il, 1]
                LabDemand       = endo.𝐨[iz, ia] ? endo.𝐥[iz, ia] : 0.0
                JobsDestroyed   = max(0.0, l⃗[il] - LabDemand)
                Switchers       = (l⃗[il] > 0) && (endo.𝐨[iz, ia] == false)

                # C2. Start adding up 
                endo.S          += Mass * Switchers
                endo.D          += Mass * JobsDestroyed
            end
        end
    end
    endo.JD     = endo.S + endo.D
end

# 2. Update labour market (unemployment and matching)
function fnUpdateLabourMarket!(params, endo)

    # A. Unpacking business 
    @unpack γ = params

    # B. Compute the key values
    endo.U = sum(@views(endo.g[:, :, :, 2]))
    endo.M = γ * (endo.U + endo.JD)
    endo.W = sum(@views(endo.g[:, :, 1, 1]))
    endo.E = 1.0 - endo.U - endo.W
end

# 3. Compute flows 
function fnComputeFlows(iz, ia, il, iu, params, endo)

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
    JDR = min(1.0, max(0.0, endo.D / endo.W))
    WU  = 0.0 + (il == 1) * JDR * (1 - JFR) * (iu == 1) * (endo.𝐨[iz, ia] == false)
    WE  = 0.0 + (il == 1) * (endo.𝐨[iz, ia] == true) * (iu == 1)
    EU  = 0.0 + (il > 1) * (endo.𝐨[iz, ia] == false) * (iu == 1) * (1 - JFR)
    EW  = 0.0 + (il > 1) * (endo.𝐨[iz, ia] == false) * (iu == 1) * JFR
    UE  = 0.0 + (endo.𝐨[iz, ia] == true) * (iu == 2)
    UW  = 0.0 + (endo.𝐨[iz, ia] == false) * (iu == 2) * JFR
    EE  = 0.0 + (il > 1) * (endo.𝐨[iz, ia] == true) * (iu == 1)
    WW  = 0.0 + (il == 1) * (endo.𝐨[iz, ia] == false) * (iu == 1) * ((1 - JDR) + (JDR * JFR))
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
    @unpack z⃗, a⃗, l⃗, δᵈⁱˢᵗ, ψ, μ⃗, Nᶻ, Nᵃ, Nˡ, Nᵘ = params

    # B. Starting the loop for a PDF
    εᵈⁱˢᵗ   = 1.0
    gⁿᵉˣᵗ   = zeros(Nᶻ, Nᵃ, Nˡ, Nᵘ)
    𝓃ᵈⁱˢᵗ   = 1
    B       = zeros(Nᵃ, Nˡ, 2)

    # C. Precompute the invariant elements 
    # I. Initialisation 
    ibᵃ     = zeros(Int, Nᶻ, Nᵃ)
    iuᵃ     = zeros(Int, Nᶻ, Nᵃ)
    wᵃ      = zeros(Float64, Nᶻ, Nᵃ)
    ibˡ     = zeros(Int, Nᶻ, Nᵃ)
    iuˡ     = zeros(Int, Nᶻ, Nᵃ)
    wˡ      = zeros(Float64, Nᶻ, Nᵃ)

    # II. Loop 
    @inbounds for ia in eachindex(a⃗)
        for iz in eachindex(z⃗)

            # (a) Next asset mass 
            Aₜ          = endo.𝐚[iz, ia]
            ibᵃ[iz,ia]  = clamp(searchsortedlast(a⃗, Aₜ), 1, length(a⃗) - 1)
            iuᵃ[iz,ia]  = ibᵃ[iz,ia] + 1
            wᵃ[iz,ia]   = clamp((a⃗[iuᵃ[iz,ia]] - Aₜ) / (a⃗[iuᵃ[iz,ia]] - a⃗[ibᵃ[iz,ia]]), 0.0, 1.0)
            # (b) Labour mass 
            Lₜ          = endo.𝐨[iz, ia] ? endo.𝐥[iz, ia] : 0.0
            ibˡ[iz,ia]  = clamp(searchsortedlast(l⃗, Lₜ), 1, length(l⃗) - 1)
            iuˡ[iz,ia]  = ibˡ[iz,ia] + 1
            wˡ[iz,ia]   = clamp((l⃗[iuˡ[iz,ia]] - Lₜ) / (l⃗[iuˡ[iz,ia]] - l⃗[ibˡ[iz,ia]]), 0.0, 1.0)
        end 
    end 

    # D. Begin the loop 
    while (εᵈⁱˢᵗ > δᵈⁱˢᵗ)

        # E. Aggregate states 
        fnJobDestruction!(params, endo)         # Getting JD, S 
        fnUpdateLabourMarket!(params, endo)     # Getting U, M, W, E  

        @inbounds for iu in 1:2 
            for il in eachindex(l⃗) 
                for ia in eachindex(a⃗)
                    for iz in eachindex(z⃗)
                        # I. Common terms
                        Mass = endo.g[iz, ia, il, iu]
                        if Mass == 0.0
                            continue
                        end

                        # II. Update the flows 
                        # (c) Get flows that matter
                        f¹², f²¹, f²², f¹¹ = fnComputeFlows(iz, ia, il, iu, params, endo)
                        FlowU = f¹² + f²²
                        FlowO = f²¹ + f¹¹

                        # III. With switching productivity
                        # To unemployment (iu=2)
                        B[ibᵃ[iz,ia], ibˡ[iz,ia], 2]    += Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowU
                        B[ibᵃ[iz,ia], iuˡ[iz,ia], 2]    += Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowU
                        B[iuᵃ[iz,ia], ibˡ[iz,ia], 2]    += Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowU
                        B[iuᵃ[iz,ia], iuˡ[iz,ia], 2]    += Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowU
                        # To employment/entrepreneurship (iu=1)
                        B[ibᵃ[iz,ia], ibˡ[iz,ia], 1]    += Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowO
                        B[ibᵃ[iz,ia], iuˡ[iz,ia], 1]    += Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowO
                        B[iuᵃ[iz,ia], ibˡ[iz,ia], 1]    += Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowO
                        B[iuᵃ[iz,ia], iuˡ[iz,ia], 1]    += Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowO

                        # IV. The same productivity
                        # (c) Update (asset, labour) = (BB, BU, UB, UU)
                        # To unemployment 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], ibˡ[iz,ia], 2]    += ψ * Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowU
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], iuˡ[iz,ia], 2]    += ψ * Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], ibˡ[iz,ia], 2]    += ψ * Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], iuˡ[iz,ia], 2]    += ψ * Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowU
                        # To others 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], ibˡ[iz,ia], 1]    += ψ * Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowO
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], iuˡ[iz,ia], 1]    += ψ * Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], ibˡ[iz,ia], 1]    += ψ * Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], iuˡ[iz,ia], 1]    += ψ * Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowO
                    end
                end
            end
        end

        # IV. Add the mass to the updated productivity states
        @inbounds for iu in 1:2, il in 1:Nˡ ,ia in 1:Nᵃ, iz in 1:Nᶻ
            gⁿᵉˣᵗ[iz, ia, il, iu] += (1 - ψ) * μ⃗[iz] * B[ia, il, iu]
        end
        
        # V. Updating business 
        εᵈⁱˢᵗ = maximum(abs(gⁿᵉˣᵗ[i] - endo.g[i]) for i in eachindex(gⁿᵉˣᵗ))
        endo.g .= gⁿᵉˣᵗ
        fill!(gⁿᵉˣᵗ, 0.0)
        fill!(B, 0.0)
        # if 𝓃ᵈⁱˢᵗ % 100 == 0
        #     println("→→→ Dist. iteration: $𝓃ᵈⁱˢᵗ, εᵈⁱˢᵗ: $(round(εᵈⁱˢᵗ, digits=8))")
        # end
        𝓃ᵈⁱˢᵗ += 1
    end
end

# 5. Aggregate states 
function fnAggregateStates!(params, endo)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, δᵈⁱˢᵗ, a⃗ = params

    # B. Iterate to find the distribution and its marginal 
    fnForwardIteration!(params, endo)
    endo.g̃ .= dropdims(sum(endo.g, dims=4), dims=4)

    # D. Update the labour market 
    fnJobDestruction!(params, endo)
    fnUpdateLabourMarket!(params, endo)

    # E. Weight for the marginal z × a distribution
    ωᵃ = dropdims(sum(endo.g, dims=(3, 4)), dims=(3, 4))

    # F. Capital demand and supply 
    endo.Kˢ = sum(ωᵃ .* a⃗')
    endo.Kᵈ = sum(ωᵃ .* endo.𝐤 .* endo.𝐨)

    # G. Labour demand and supply 
    endo.Lˢ = endo.W
    endo.Lᵈ = sum(ωᵃ .* endo.𝐥 .* endo.𝐨)
end

# -----------------------------------------------
# C. Solving  
# -----------------------------------------------

# 0. Convex updating 
function fnConvexUpdating(f, guesses::Tuple, method=nothing; 
                          loading = 0.2, xatol=1e-5, max_iter=250, 
                          init=nothing, kwargs...)
    
    # A. Start 
    xˣ          = (isnothing(init) || init == 0.0) ? (guesses[1] + guesses[2]) / 2.0 : float(init)
    prev_res    = 0.0
    step        = loading 

    for iter in 1:max_iter
        residual= f(xˣ)
        
        # B. Break immediately when tolerance is met!
        if abs(residual) <= xatol
            return xˣ
        end

        # C. Minor acceleration if needed 
        if iter > 1
            if sign(residual) != sign(prev_res)
                step *= 0.5
            else
                step = min(step * 1.1, 0.1) 
            end
        end

        # D. Update and clamp to bounds
        x_new       = xˣ + (residual * step)
        xˣ          = clamp(x_new, guesses[1], guesses[2])
        prev_res    = residual
    end

    @warn "Convex updating hit the ceiling ($max_iter). Final error: $(abs(prev_res))"
    return xˣ
end

# 1. Labour residual.
function fnLabourResidual!(w, params, endo)
    # A. Update the state
    endo.wₜ = w

    # B. Run the parallelised engine
    fnVFI!(params, endo)
    fnAggregateStates!(params, endo)
    
    # C. Return the labour market error
    εᴸ      = (endo.Lᵈ / max(endo.Lˢ, 1e-4)) - 1.0
    println("→→→ Wage search        | w = $(round(w, digits=4)), εᴸ = $(round(εᴸ, digits=7))")
    return εᴸ
end

# D. A struct to hold the variables needed for the labor residual
struct LabourResidualObjective{P, E}
    params::P
    endo::E
end

# E. Make it callable
function (obj::LabourResidualObjective)(w)
    return fnLabourResidual!(w, obj.params, obj.endo)
end

# 3. Find the clearing tax
function fnClearingTax!(params, endo)

    # A. Unpacking business 
    @unpack w̲, w̅, δᴸ,τ̲, τ̅, δᵗ,ηᵗ = params

    # B. Prepare the loop 
    εᵗ      = 10.0
    
    # C. Start the loop 
    while abs(εᵗ) > δᵗ

        # C1. Get the wage
        Oᴸ          = LabourResidualObjective(params, endo)
        ϵ           = Oᴸ(endo.wₜ)
        wˣ = try
        if abs(ϵ) <= δᴸ
            endo.wₜ
        elseif ϵ > 0
            find_zero(Oᴸ, (endo.wₜ, w̅), Bisection(), xatol=δᴸ)
        else
            find_zero(Oᴸ, (w̲, endo.wₜ), Bisection(), xatol=δᴸ)
        end
        catch
            find_zero(Oᴸ, (w̲, w̅), Bisection(), xatol=δᴸ)
        end
        endo.wₜ     = wˣ

        # C2. Error & update
        εᵗ          = endo.wₜ * endo.U - endo.τₜ
        endo.τₜ     = ηᵗ * endo.wₜ * endo.U + (1 - ηᵗ) * endo.τₜ
    
        # D. Warm message :) 
        println("→→ Tax loop            | τ = $(round(endo.τₜ, digits=4)), εᵗ = $(round(εᵗ, digits=7)) [Cleared w = $(round(endo.wₜ, digits=4))]")
    end
end 

# 4. Capital residual 
function fnCapitalResidual!(params, endo, r)
        
        # A. Unpacking business 
        @unpack w̲, w̅, δᴸ = params

        # B. Solve 
        endo.rₜ = r 
        fnClearingTax!(params, endo)

        # C. Return error
        εᴷ      = (endo.Kᵈ / max(endo.Kˢ, 1e-4)) - 1.0
        println("⋆ Capital loop     | r = $(round(r, digits=4)), εᴷ = $(round(εᴷ, digits=4)) [Cleared τ = $(round(endo.τₜ, digits=4))]")
        return εᴷ
    end

# D. Struct for the capital residual
struct CapitalResidualObjective{P, E}
    params::P
    endo::E
end

# E. Make it callable
function (obj::CapitalResidualObjective)(r)
    return fnCapitalResidual!(obj.params, obj.endo, r)
end

# 5. Solve for steady state 
function fnSolveSteadyState!(params, endo)
    
    # A. Unpacking business 
    @unpack r̲, r̅, δʳ = params

    # B. The final solve
    if endo.wₜ == 0.0; endo.wₜ = 1.532166; end
    if endo.τₜ == 0.0; endo.τₜ = 0.101028; end
    if endo.rₜ == 0.0; endo.rₜ = 0.009158; end
    Oᶜ          = CapitalResidualObjective(params, endo)
    rˣ          = fnConvexUpdating(Oᶜ, (r̲, r̅), xatol = δʳ, init = endo.rₜ)
    endo.rₜ     = rˣ

    # C. Lock in the results
    fnVFI!(params, endo)
    fnAggregateStates!(params, endo)
    println("\n--- Steady state ---")
    println("Wage (w):   $(round(endo.wₜ, digits=6))")
    println("Interest(r):$(round(endo.rₜ, digits=6))")
    println("Tax (τ):    $(round(endo.τₜ, digits=6))")
end
 
# 6. Print the results 

function fnPrintCalibrationElements(params, endo)
    # A. Unpack parameters to keep code clean
    @unpack Nᶻ, Nᵃ, a⃗, ψ = params
    
    # B. Obtain the marginal distribution over (z, a)
    ωᵃ = dropdims(sum(endo.g, dims=(3, 4)), dims=(3, 4))
    
    # C. Standard targets 
    model_r         = endo.rₜ
    t_e             = sum(ωᵃ .* endo.𝐨)
    model_exit      = endo.S / max(t_e, 1e-8)
    ext_finance     = sum(ωᵃ .* max.(endo.𝐤 .- a⃗', 0.0) .* endo.𝐨)
    model_ext_fin   = ext_finance / max(endo.Kᵈ, 1e-8)
    
    # D. Employment statistics 
    emp_levels      = endo.𝐥[endo.𝐨]       # Labor demand for those firms
    emp_masses      = ωᵃ[endo.𝐨]           # Mass of those firms
    sort_idx_emp    = sortperm(emp_levels, rev=true)
    sorted_emp      = emp_levels[sort_idx_emp]
    sorted_mass_emp = emp_masses[sort_idx_emp]
    target_ent_mass = 0.10 * sum(sorted_mass_emp)
    total_emp       = sum(sorted_emp .* sorted_mass_emp)
    top10_emp_sum   = 0.0
    cum_mass_emp    = 0.0
    for i in eachindex(sorted_mass_emp)
        if cum_mass_emp + sorted_mass_emp[i] <= target_ent_mass
            top10_emp_sum   += sorted_emp[i] * sorted_mass_emp[i]
            cum_mass_emp    += sorted_mass_emp[i]
        else
            # Grab the fractional piece that crosses the 10% boundary
            rem_mass        = target_ent_mass - cum_mass_emp
            top10_emp_sum   += sorted_emp[i] * rem_mass
            break
        end
    end
    model_top10_emp         = total_emp > 0.0 ? (top10_emp_sum / total_emp) : 0.0

    # E. Top 5% earnings share
    earnings_levels         = zeros(Nᶻ * Nᵃ)
    earnings_masses         = zeros(Nᶻ * Nᵃ)    
    idx = 1
    for ia in 1:Nᵃ, iz in 1:Nᶻ
        earnings_levels[idx] = endo.𝐨[iz, ia] ? endo.Π[iz, ia] : endo.wₜ
        earnings_masses[idx] = ωᵃ[iz, ia]
        idx += 1
    end
    sort_idx_earn    = sortperm(earnings_levels, rev=true)
    sorted_earn      = earnings_levels[sort_idx_earn]
    sorted_mass_earn = earnings_masses[sort_idx_earn]
    target_pop_mass  = 0.05 * sum(ωᵃ) 
    total_earn       = sum(sorted_earn .* sorted_mass_earn)
    top5_earn_sum    = 0.0
    cum_mass_earn    = 0.0
    for i in eachindex(sorted_mass_earn)
        if cum_mass_earn + sorted_mass_earn[i] <= target_pop_mass
            top5_earn_sum   += sorted_earn[i] * sorted_mass_earn[i]
            cum_mass_earn   += sorted_mass_earn[i]
        else
            rem_mass        = target_pop_mass - cum_mass_earn
            top5_earn_sum   += sorted_earn[i] * rem_mass
            break
        end
    end
    model_top5_earn         = total_earn > 0.0 ? (top5_earn_sum / total_earn) : 0.0

    # F. Print results 
    println("\nTable 1\nCalibration.")
    println(repeat("-", 90))
    @printf("%-50s %-10s %-16s %-10s\n", "", "US data", "Original model", "Current model")
    println(repeat("-", 90))
    @printf("%-50s %-10.2f %-16.2f %-10.2f\n", "Top 10% employment", 0.69, 0.69, model_top10_emp)
    @printf("%-50s %-10.2f %-16.2f %-10.2f\n", "Top 5% earnings share", 0.30, 0.30, model_top5_earn)
    @printf("%-50s %-10.2f %-16.2f %-10.2f\n", "Establishment exit rate (annual)", 0.10, 0.10, model_exit)
    @printf("%-50s %-10.2f %-16.2f %-10.4f\n", "Real interest rate (annual)", 0.02, 0.02, model_r)
    @printf("%-50s %-10.2f %-16.2f %-10.4f\n", "Credit market instruments to non-financial assets", 0.70, 0.70, model_ext_fin)
    println(repeat("-", 90))
end

