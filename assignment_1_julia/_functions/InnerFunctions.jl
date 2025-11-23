##          Inner functions

#           Content:
#           1.      Intensive labour supply
#           2A.     Working value 
#           2B.     Non-working value
#           3.      Extensive labour supply
#           4.      Aggregate labour supply
#           5.      Aggregate labour demand
#           6.      Government lump sum


##################################################
##          1. Intensive labour supply          ##
##################################################

function    fnIntensiveLabourSupply(w,r,T,z,a,τ,η,χ,β)

    #       Function settings
    iConsumption    = 0.6751
    iConsumptionNext=0.0
    iConvergence    = false
    iIteration      = 1
    pErrorTol       = 1e-5
    pStepSize       = 0.9

    #       Start the loop
    while iIteration < 1000

        #   From Euler
        iRHS                = z * w * (1-τ) * (z * w * (1-τ) / (η * iConsumption))^χ + a * (1 + r * (1 - τ)) + T

        #   Adjust if needed
        if iRHS     < 0.0
            iConsumption    = 0.5 * iConsumption
            continue
        end 

        #   Update consumption & compute the error
        iConsumptionNext    = iRHS / (1 + β) 
        iError              = abs(iConsumptionNext-iConsumption)

        #   Convergence criterion 
        if iError < pErrorTol
            iConvergence    = true
            break
        else 
            # Update consumption
            iConsumption    = pStepSize * iConsumption + (1 - pStepSize) * iConsumptionNext
        end 
        iIteration          = iIteration + 1
    end 

    #       Erorr for a failed loop
    if !iConvergence
        println("w = ", w, " z = ",z," T = ", T, " η = ",η," χ = ",χ, " β = ",β)
        throw("Warnning: Convergence procedure has failed.")
    end  
    
    #       Save n
    n           =  ((z * w * (1-τ)) / (iConsumptionNext * η))^χ
    return iConsumption, n
end

##################################################
##          2A. Working value function          ##
##################################################

function    fnWorkingValue(c,n,β,χ,η)

    #       Assets & Frisch elasticity
    ap      = β * c
    Frisch  = 1.0 + (1.0 / χ)

    #       Results 
    W       = log(c) - η * (1 / Frisch) * n^(Frisch) + β * log(ap)
    return W
end

##################################################
##          2B. Non-working value function      ##
##################################################

function    fnNonWorkingValue(T,b,a,r,τ,β)

    #       Consumption from the non-working Euler
    Budget  = b + a * (1.0+r * (1 - τ)) + T 
    c       = Budget  / (1.0 + β)

    #       Assets 
    ap      = β * c 

    #       Value function
    NW      = log(c) + β * log(ap)
    return NW
end 

##################################################
##          3. Extensive labour supply          ##
##################################################

function    fnExtensiveLabourSupply(w,r,T,z,a,b,τ,η,χ,β)

    #       Compute intensive labour supply for working HHs
    cₑ, nₑ  = fnIntensiveLabourSupply(w,r,T,z,a,τ,η,χ,β)

    #       Value functions
    W       = fnWorkingValue(cₑ,nₑ,β,χ,η)
    N       = fnNonWorkingValue(T,b,a,r,τ,β)

    #       Choose working
    if W >= N 
        Decision    = 1.0
        Hours       = nₑ
    else
        Decision    = 0.0
        Hours       = 0.0
    end

    #       Produce results 
    return Decision, Hours
end 

##################################################
##          4. Aggregate labour supply          ##
##################################################

function    fnAggregateLabourSupply(w,r,T,z̄,a,b,τ,η,χ,β,σ)

    #       Set the distribution
    DistributionZ   = LogNormal(log(z̄) - 0.5 * σ^2, σ)

    #       Labour supply 
    LabourSupply    = solve(IntegralProblem((x,p) ->  pdf(DistributionZ,x) * x * fnExtensiveLabourSupply(w,r,T,x,a,b,τ,η,χ,β)[2],(0.0,Inf)),QuadGKJL()).u
    Employment      = solve(IntegralProblem((x,p) ->  pdf(DistributionZ,x) * fnExtensiveLabourSupply(w,r,T,x,a,b,τ,η,χ,β)[1],(0.0,Inf)),QuadGKJL()).u
    LabourSupply2M  = solve(IntegralProblem((x,p) ->  pdf(DistributionZ,x) * (x * fnExtensiveLabourSupply(w,r,T,x,a,b,τ,η,χ,β)[2])^2,(0.0,Inf)),QuadGKJL()).u

    #       Print results
    return LabourSupply, Employment, LabourSupply2M
end

##################################################
##          5. Aggregate labour demand          ##
##################################################

function    fnAggregateLabourDemand(K,L,A,α)
    w       = A * (1.0 - α) * (K / L)^α
    return w
end

##################################################
##          6. Government lump transfer         ##
##################################################

function    fnGovernmentLumpTransfer(L,E,w,a,r,τ,b)
    T       = (L * w + a * r) * τ - (1.0 - E) * b 
    return T
end