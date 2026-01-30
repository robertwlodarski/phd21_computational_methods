# Model infrastructure

# Content
# 1. Endogenous labour grid 
# 2. Bargained wage
# 3. Initial value function guest  
# 4. Value function iteration 
# 5. Aggregation 

# 1. Endogenous labour grid 

function fn⃗(params::ModelParameters,p,f,q)
    # A. Unpacking business 
    @unpack η, α, x̲, b, β, c, x̅, s̅ₙ, s̲ₙ, Δₙ,N₁,N₂,N₃,N₄= params

    # B. Define n̲ and n̅
    n̲       =  s̲ₙ * ( ((1 - η*(1-α)) / (p*x̲*α)) * (b + (1/(1-η)) * (η*β*c*(f/q) + (1 - η*(1-α))*c/q)) )^(-1/(1-α))
    n̅       = s̅ₙ * ( (1 - η * (1 - α)) / (p * x̅ * α) * (b + (1 / (1 - η)) *  (η * β * f * c / q)) )^(- 1 / (1 - α))

    # C. Define the n grid 
    return exp.(collect(log(n̲) : Δₙ : log(n̅)))     
end 

# 2. Bargained wage

function fW(params::ModelParameters,p,f,q,n⃗)
    # A. Unpacking business 
    @unpack η, α, b, β, c, x⃗ = params 

    # B. Compute wage 
    W    = η * ( p .* x⃗ .* α .* n⃗'.^(α - 1) ./ (1 - η * (1 - α)) + β * f * c / q) + (1 - η) * b 

    # C. Return 
    return W
end 

# 3. Initial value function guess 

function fΠ⁰!(params::ModelParameters,endo::EndogenousVariables,p,n⃗,W,q)
    # A. Unpacking business 
    @unpack α, Nₓ, x⃗, πˢᶜᵃˡᵉ, β, c, λ, W⃗ₓ = params 

    # B. Compute flow profit & initial guess 
    endo.Πᶠˡᵒʷ      = p .* x⃗ .* (n⃗' .^(α)) .- W .*  n⃗'
    endo.Πᶜ         = πˢᶜᵃˡᵉ * (1 / (1 - β)) .* Πᶠˡᵒʷ

    # C. Loop for every x 
    vⁿᵉʷ            = zeros(Nₓ)
    nˣ              = zeros(Nₓ)
    for i in 1:Nₓ
        val, id     = findmax(view(Πᶜ,i,:))
        vⁿᵉʷ[i]     = val 
        nˣ[i]       = n⃗[id]
    end 

    # D. Firing and hiring functions 
    Πᶠ              = vⁿᵉʷ .* (nˣ < n⃗') - 1e8 * (nˣ > n⃗')
    Πʰ              = (vⁿᵉʷ .- c / q .* (nˣ - n⃗')) .* (nˣ > n⃗') - 1e8 * (nˣ < n⃗')

    # E. Update the key value functions 
    endo.Π          .= max.(Πᶠ,max.(Πʰ,Πᶜ))
    endo.𝔼Π         .= (1 - λ) .* endo.Π + λ .* W⃗ₓ' * endo.Π 
end 

# 4. Value function iteration 
function fVFI!(params::ModelParameters,endo::EndogenousVariables,p,f,q)

    # A. Unpacking business
    @unpack Nₓ, β, c, λ, W⃗ₓ, n̅ˢ,N₁,N₂,N₃,N₄, x⃗, x̲, N̅₁, N̅₂, x̅ = params

    # B. Construct the employment grid, n⃗, and matrix of wages 
    n⃗       = fn⃗(params,p,f,q)
    W       = fW(params,p,f,q,n⃗)

    # C. Initial value function guesses
    fΠ⁰!(params,endo,p,n⃗,W,q)
    Πᵒˡᵈ    = endo.Π
    𝔼Π      = endo.𝔼Π
    Πᶠˡᵒʷ   = endo.Πᶠˡᵒʷ

    # D. Initialise values and start the loop
    Πᶠ      = zeros(size(Πᵒˡᵈ)) 
    Πʰ      = zeros(size(Πᵒˡᵈ))
    Πⁿᵉʷ    = zeros(size(Πᵒˡᵈ))
    vᶠ      = zeros(Nₓ)
    nᶠ      = zeros(Nₓ)
    vʰ      = zeros(Nₓ)
    nʰ      = zeros(Nₓ)
    𝕀ˢᵗᵒᵖ   = false 
    𝕀ᶠᵃˢᵗ   = true 
    nᵛ      = 1
    nˢ      = 1
    while 𝕀ˢᵗᵒᵖ == false 
        Πᶜ  = Πᶠˡᵒʷ .+ β .* 𝔼Π

        # E. Coarse updating 
        if 𝕀ᶠᵃˢᵗ   == true 

                # i. Compute the values of firing and hiring
                for i in 1:Nₓ
                    # Fire 
                    valᶠ, idᶠ   = findmax(view(Πᶜ,i,:))
                    vᶠ[i]       = valᶠ 
                    nᶠ[i]       = n⃗[idᶠ]
                    # Hire 
                    valʰ, idʰ   = findmax(view(Πᶜ,i,:).-(c/q).*n⃗)
                    vʰ[i]       = valʰ 
                    nʰ[i]       = n⃗[idʰ]
                end 
                Πᶠ              .= vᶠ .* (nᶠ < n⃗') - 1e8 * (nᶠ > n⃗')
                Πʰ              .= (vʰ .+ c / q .* n⃗') .* (nʰ > n⃗') - 1e8 * (nʰ < n⃗')

                
                # ii. Update values and error terms 
                Πⁿᵉʷ            .= max.(Πᶠ,max.(Πʰ,Πᶜ))
                𝔼Π              .= (1 - λ) .* Πⁿᵉʷ + λ .* W⃗ₓ' * Πⁿᵉʷ 
                εᵛᶠⁱ            = maximum(abs.((Πⁿᵉʷ-Πᵒˡᵈ)./Πᵒˡᵈ))
                nᵛ              +=1
                Πᵒˡᵈ            .= Πⁿᵉʷ
        end 

        # F. Refined updating
        if 𝕀ᶠᵃˢᵗ   == false 
            
            # i. Compute values of firing and hiring 
            for i in 1:Nₓ
                # Fire 
                ℑᶠ              = CubicSplineInterpolation(n⃗,view(Πᶜ,i,:))
                ℜᶠ              = optimize(n -> -ℑᶠ(n),n⃗[1],n⃗[end])
                nᶠ[i]           = Optim.minimizer(ℜᶠ)
                vᶠ[i]           = -Optim.minimum(ℜᶠ)
                # Hire 
                ℑʰ              = CubicSplineInterpolation(n⃗,view(Πᶜ,i,:))
                ℜʰ              = optimize(n -> -ℑʰ(n)+(c/q)*n,n⃗[1],n⃗[end])
                nʰ[i]           = Optim.minimizer(ℜʰ)
                vʰ[i]           = -Optim.minimum(ℜʰ)
            end 
            # Compute 
            Πᶠ                  .= vᶠ .* (nᶠ < n⃗') - 1e8 * (nᶠ > n⃗')
            Πʰ                  .= (vʰ .+ c / q .* n⃗') .* (nʰ > n⃗') - 1e8 * (nʰ < n⃗')
            
            # ii. Update values and error terms 
            Πⁿᵉʷ                .= max.(Πᶠ,max.(Πʰ,Πᶜ))
            𝔼Π                  .= (1 - λ) .* Πⁿᵉʷ + λ .* W⃗ₓ' * Πⁿᵉʷ 
            εᵛᶠⁱ                = maximum(abs.((Πⁿᵉʷ-Πᵒˡᵈ)./Πᵒˡᵈ))
            # Stop when too many splines 
            if nˢ == n̅ˢ  
                𝕀ˢᵗᵒᵖ = true 
            end 
            nᵛ                  +=1
            nˢ                  +=1
            Πᵒˡᵈ                .= Πⁿᵉʷ
        end 

        # G. Grid refinement 
        if (εᵛᶠⁱ < δʳᵉᶠ && 𝕀ᶠᵃˢᵗ == true)
            # i. Save the old value  
            n⃗ᵒˡᵈ                = copy(n⃗) 

            # ii. Define grid boundaries 
            n̲₁                  = nʰ[1]
            n̲₂                  = 0.975*nᶠ[1]
            n̅₁                  = 1.025*nʰ[Nₓ]
            n̅₂                  = nᶠ[Nₓ]

            # ii. Grids 
            ñ⃗₁                  = 10 .^ range(log10(n̲₁),log10(n̲₂),length=N₁)
            ñ⃗₂                  = 10 .^ range(log10(n̲₂),log10(1.25*n̲₂),length=N₂)
            ñ⃗₃                  = 10 .^ range(log10(1.25*n̲₂),log10(n̅₁),length=N₃)
            ñ⃗₄                  = 10 .^ range(log10(n̅₁),log10(n̅₂),length=N₃)
            n⃗                   = unique([ñ⃗₁;ñ⃗₂;ñ⃗₃;ñ⃗₄])
            Nₙ                  = length(n⃗)

            # iv. Interpolation station
            W                   = fW(params,p,f,q,n⃗)
            Πᶠˡᵒʷ               = p .* x⃗ .* (n⃗' .^(α)) .- W .*  n⃗'
            𝔼Πⁿᵉʷ               = zeros(Nₓ,length(n⃗))
            for i in 1:Nₓ
                ℑʳ              = CubicSplineInterpolation(n⃗ᵒˡᵈ,view(𝔼Π,i,:))
                𝔼Πⁿᵉʷ[i,:]      = ℑʳ(n⃗)
            end
            𝔼Π                  = 𝔼Πⁿᵉʷ

            # v. Settings 
            𝕀ᶠᵃˢᵗ               = false 
            nˢ                  += 1
            Πᵒˡᵈ                = copy(Πᶠˡᵒʷ) 
            Πᶠ                  = zeros(Nₓ, Nₙ)
            Πʰ                  = zeros(Nₓ, Nₙ)
            Πⁿᵉʷ                = zeros(Nₓ, Nₙ)
        end 
    end  

    # 5. Produce the policy functions of interest
    # A. Employment policy
    ℑⁿˡ     = CubicSplineInterpolation(x⃗,nᶠ)
    n̲ᵖ      = ℑⁿˡ(1.001 * x̲)
    n̅ᵖ      = nʰ[Nₓ]
    n̂ᵖ      = 1.25 * n̲ᵖ
    n̂⃗₁      = 10 .^ range(log10(n̲ᵖ),log10(n̂ᵖ),length=N̅₁)
    n̂⃗₂      = 10 .^ range(log10(n̂ᵖ),log10(n̅ᵖ),length=N̅₂)
    endo.n⃗  = unique([n̂⃗₁;n̂⃗₂])

    # B. Indices for firing and hiring thresholds 
    𝓃₁      = findlast(nᶠ .< n̅ᵖ)
    𝓃₂      = findfirst(nʰ .> n̲ᵖ)

    # C. Firing threshold 
    𝕟ᴿ      = nᶠ[1:𝓃₁]
    𝒾ᴿ      = unique(i -> 𝕟ᴿ[i],1:length(𝕟ᴿ)) 
    𝕩ᴿ      = x⃗[1:𝓃₁]
    ℑᴿ      = CubicSplineInterpolation(𝕟ᴿ[𝒾ᴿ],𝕩ᴿ[𝒾ᴿ])
    endo.R⃗  = ℑᴿ(endo.n⃗)
    endo.∂R⃗ = Interpolations.derivative.(Ref(ℑᴿ),endo.n⃗)

    # D. Hiring threshold 
    𝕟ᴿⱽ     = nʰ[𝓃₂:end]
    𝒾ᴿⱽ     = unique(i -> 𝕟ᴿⱽ[i],1:length(𝕟ᴿⱽ))
    𝕩ᴿⱽ     = x⃗[𝓃₂:end]
    ℑᴿⱽ     = CubicSplineInterpolation(𝕟ᴿⱽ[𝒾ᴿⱽ],𝕩ᴿⱽ[𝒾ᴿⱽ])
    endo.R⃗ᵥ = min.(ℑᴿⱽ(endo.n⃗),x̅)
    endo.∂R⃗ᵥ= Interpolations.derivative.(Ref(ℑᴿⱽ),endo.n⃗)
    
end 

# 5. Aggregation 
function fAggregation!(params::ModelParameters,endo::EndogenousVariables,p,f,q)

    # A. Unpacking business 
    @unpack x̲, x⃗, ξ, p̄ₓ = params 

    # B. Compute CDFs, PDFs, and expectation 
    𝐆R⃗ᵥ     = (1 .- (x̲ ./ endo.R⃗ᵥ).^ξ) ./ p̄ₓ
    𝐠R⃗ᵥ     = ((1 / p̄ₓ) * ξ * x̲^ξ) ./ ((endo.R⃗ᵥ).^(ξ+1)) 
    𝐆R⃗      = (1 .- (x̲ ./ endo.R⃗).^ξ) ./ p̄ₓ
    𝐠R⃗      = ((1 / p̄ₓ) * ξ * x̲^ξ) ./ ((endo.R⃗).^(ξ+1)) 
    𝐇n⃗      = 𝐆R⃗ ./ (1 - 𝐆R⃗ᵥ + 𝐆R⃗)
    𝐡n⃗      = ((1 - 𝐆R⃗ᵥ) .* 𝐠R⃗ .* endo.∂R⃗ + 𝐆R⃗ .* 𝐠R⃗ᵥ .* endo.∂R⃗ᵥ) ./ ((1 - 𝐆R⃗ᵥ + 𝐆R⃗).^2)
    𝔼x      = x̲^ξ * (ξ /(ξ - 1)) * (endo.R⃗.^(-ξ+1)-endo.R⃗ᵥ.^(-ξ+1)) ./ (𝐆R⃗ᵥ .- 𝐆R⃗)
    
end 