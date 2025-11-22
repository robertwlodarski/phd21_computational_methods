# Main functions

##################################################
##          1. Wage solver                      ##
##################################################

function        fnWageSolver(A,α,r,z̄,a,b,τ,η,χ,β,σ)

    #       Initialise 
    iWage       = 0.99
    iT          = 0.04
    iWageNext   = 0
    iTNext      = 0
    iError      = 1
    iIteration  = 0
    iConvergence= false
    pMaxIter    = 1000
    pStepSize   = 0.7

    #       Starting the process
    while iIteration < pMaxIter
        
        #   Labour supply
        L, E, _         = fnAggregateLabourSupply(iWage,r,iT,z̄,a,b,τ,η,χ,β,σ)

        #   Labour demand (and transfers)
        iWageNext       = fnAggregateLabourDemand(a,L,A,α)
        iTNext          = fnGovernmentLumpTransfer(L,E,iWage,a,r,τ,b)

        #   Update error 
        iError          = abs(iWageNext - iWage) + abs(iTNext - iT)

        #   Convergence station
        if iError < 1e-5
            iConvergence= true
            break
        else
            iWage       = pStepSize * iWage + (1 - pStepSize) * iWageNext
            iT          = pStepSize * iT    + (1 - pStepSize) * iTNext  
        end 
        iIteration      = iIteration + 1
    #       End the while loop
    end 

    #       Print messages
    println("Solver has finished successfully.")
    println("Error              = ", iError)
    println("# of iterations    = ", iIteration)
    println("Wage               = ", iWageNext)
    println("Lump trasnfer      = ", iTNext)

    #       Save results
    return iWageNext, iTNext
end 