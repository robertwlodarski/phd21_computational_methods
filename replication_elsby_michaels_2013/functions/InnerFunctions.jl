# Model infrastructure

# Content
# 1. Endogenous labour grid 
# 2. Bargained wage
# 3. Flow profit 
# 4. Initial profit guess 

# 1. Endogenous labour grid 

function fnEndogenousLabourGrid(params::ModelParameters,p,f,q)
    # A. Unpacking business 
    @unpack η, α, x̲, b, β, c, x̅, s̅ₙ, s̲ₙ, Δₙ = params

    # B. Define n̲ and n̅
    n̲       =  s̲ₙ * ( ((1 - η*(1-α)) / (p*x̲*α)) * (b + (1/(1-η)) * (η*β*c*(f/q) + (1 - η*(1-α))*c/q)) )^(-1/(1-α))
    n̅       = s̅ₙ * ( (1 - η * (1 - α)) / (p * x̅ * α) * (b + (1 / (1 - η)) *  (η * β * f * c / q)) )^(- 1 / (1 - α))

    # C. Define the n grid 
    return exp.(collect(log(n̲) : Δₙ : log(n̅)))     
end 

# 2. Bargained wage

function fnBargainedWage(params::ModelParameters,p,f,q,n⃗,x⃗)
    # A. Unpacking business 
    @unpack η, α, b, β, c = params 

    # B. Compute wage 
    W    = η * ( p .* x⃗ .* α .* n⃗'.^(α - 1) ./ (1 - η * (1 - α)) + β * f * c / q) + (1 - η) * b 

    # C. Return 
    return W
end 

# 3. Flow profit 

function fnFlowProfit(params::ModelParameters,p,n⃗,x⃗,W)
    # A. Unpacking business 
    @unpack α, Nₓ = params 

    # B. Compute flow profit 
    Πᶠˡᵒʷ   = p .* x⃗ .* (n⃗' .^(α)) .- W .*  n⃗'
    return  p .* x⃗ .* (n⃗' .^(α)) .- W .*  n⃗'
end 

# 4. Initial profit guess 