# Main functions

# 1. Plain DMP w/ exogenous separation rate (naive)
# 2. Plain DMP w/ exogenous separation rate (efficient)

##################################################
##          1. Plain DMP (naive)                ##
##################################################

# Matlab style of thinking
function        fnPlainDMPSteadyStateNaive(A,p)

    ## 1.               Prepare the setting
    # A.                Unpack parameters
    @unpack μ,ξ,η,κ,b,β,λ = p

    # B.                Loop settings 
    iWeightOld          = 0.95
    iErrorGE            = 10.0
    iTolGE              = 1e-3
    θ                   = 0.0
    w                   = 0.0

    # C.                 Initial guess
    iJobFillingRate     = 0.5 

    ## 2.               Start the GE loop
    while iErrorGE > iTolGE
        
        # A.                Calculate the implied θ
        θ                   = (iJobFillingRate / μ)^(-1 / ξ)

        # B.                Calculate the implied wage
        w                   = η*(A + κ *θ) + (1 - η)*b

        # C.                Calculate the implied tightness and job finding rate 
        θnext               = ((A - w) * μ * β / (κ*(1 - β*(1- λ))))^(1 / ξ)
        iJobFillingRateNext = μ * θnext^(-ξ) 

        # D. Error term
        iErrorGE            = abs(iJobFillingRate-iJobFillingRateNext)
        iJobFillingRate     = iWeightOld * iJobFillingRate + (1 - iWeightOld) * iJobFillingRateNext
    end
    return θ, w
end 

##################################################
##          2. Plain DMP (efficient)            ##
##################################################

function        fnPlainDMPSteadyStateRoot(A,p)

    ## 1.               Prepare the setting
    # A.                Unpack parameters
    @unpack μ,ξ,η,κ,b,β,λ = p

    ## 2.               Solve the model 
    θ_sol               = find_zero(x -> fnPlainDMPResidual(x,A,p),(1e-4,50))

    ## 3.               Solve for wage 
    w_sol               = η * (A + κ * θ_sol) + (1 - η) * b    
    return θ_sol, w_sol
end 

