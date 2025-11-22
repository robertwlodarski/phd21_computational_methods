# PHD21 Computational methods: Assignment 1
# Rob Włodarski 

## 1.       Load packages
using Plots
using Distributions
using Integrals
using Optim
using NLsolve
using Roots
using LinearAlgebra
using Infiltrator

## 2.       Inner functions
#           Content:
#           1) Intensive labour supply
#           2) Extensive labour supply

function    fnIntensiveLabourSupply(w,r,T,z,a,τ,η,χ,β)
    #       Function settings
    iConsumption    = 0.6751
    iConvergence    = false
    iIteration      = 1
    pErrorTol       = 1e-5
    pStepSize       = 0.9

    #       Start the loop
    while iter < 1000

        #   From Euler
        iRHS                = z * w * (1-τ) * (z * w * (1-τ) / (η * iCons))

        #   Adjust if needed
        if iRHC     < 0
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
        println("w = ", w, "z = ",z,"T = ", T, "η = ",η)# FINISH HERE
    end  

    # FINISH HERE

    

end

## 3.       Main functions 
#           Content:
#           1) Aggregate labour supply
#           2) Labour demand (aggregate)
#           3) Government lump sum (implied)