# Model infrastructure

# Content
# 1. Endogenous labour grid 
# 2. Bargained wage
# 3. Initial value function guest  

# 1. Endogenous labour grid 

function fn⃗(params::ModelParameters,p,f,q)
    # A. Unpacking business 
    @unpack η, α, x̲, b, β, c, x̅, s̅ₙ, s̲ₙ, Δₙ = params

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
    Πᶠˡᵒʷ           = p .* x⃗ .* (n⃗' .^(α)) .- W .*  n⃗'
    Πᶜ              = πˢᶜᵃˡᵉ * (1 / (1 - β)) .* Πᶠˡᵒʷ

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

