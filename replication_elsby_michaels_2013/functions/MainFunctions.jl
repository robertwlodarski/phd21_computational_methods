# Main functions
# Explainer: Functions that solve for SS or aggregate uncertainty case 

# Content:
# 1. Compute the steady state equilibrium


# 1. Compute the steady state equilibrium
function fSteadyState!(params::ModelParameters,endo::EndogenousVariables,p)

    # A. Unpacking business 
    @unpack q̅,q̲         = params 

    # B. Define the objective function 
    r⃗           = (q̲,q̅) 
    𝒻(q)        = fEqResidual(q,p,params, endo)

    # C. Solve the model 
    println("Solving for the steady state equilibrium")
    q̂           = find_zero(𝒻,r⃗)
    println("Equilibrium job filling rate found: $q̂")

    # D. Update results 
    f̂           = fUpdatedJobFindingRate(q̂,params)
    fVFI!(params, endo,p,f̂,q̂)
    N̂,Ŷ,Ŝ,M̂,Â   = fAggregation(params, endo, p, f̂,q̂) 
    
    # E. Save the results into the structure 
    endo.Y, endo.S, endo.M, endo.A = Ŷ, Ŝ, M̂, Â
    endo.q, endo.f = q̂, f̂
    endo.U = Ŝ / f̂
    endo.N = N̂

    # E. Print the headline results  
    println("===============================================")
    println("Steady state equilibrium results:")
    println("Employment: $N̂")
    println("Production: $Ŷ")
end 


