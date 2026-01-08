##          Inner functions

#           Content:
#           1.      Next period's wealth optimsation 
#           2.      Interpolate value function


##################################################
##          1. Next wealth optimisation         ##
##################################################

function        fnWealthNext(iBudget,vMinWealth,iExpValueZ,p,g)

    # Unpacking 
    @unpack vGridA1 = g
    @unpack β       = p

    ## 1.       Define the objective function
    function    fnObjective(ap)
        ib      = searchsortedlast(vGridA1, ap)
        ib      = clamp(ib, 1, length(vGridA1)-1)
        iu      = ib + 1
        w       = (vGridA1[iu] - ap) / (vGridA1[iu] - vGridA1[ib])
        EV      = w * iExpValueZ[ib] + (1-w) * iExpValueZ[iu]
        return -(log(iBudget - ap) + β * EV)
    end 

    ## 2.       Optimise 
    iRes        = optimize(fnObjective, vMinWealth, iBudget, Brent(), abs_tol = 1e-9)
    return      Optim.minimizer(iRes) 
end

##################################################
##          2. Interpolate value function       ##
##################################################

function        fnInterpolateValueFunction(iWealthNext,iBudget,iExpValueZ,p,g)
    
    # Unpacking
    @unpack vGridA1 = g
    @unpack β       = p

    ## 1.       Get weight 
    ib      = searchsortedlast(vGridA1, iWealthNext)
    ib      = clamp(ib, 1, length(vGridA1)-1)
    iu      = ib + 1
    w       = (vGridA1[iu] - iWealthNext) / (vGridA1[iu] - vGridA1[ib])
    EV      = w * iExpValueZ[ib] + (1-w) * iExpValueZ[iu]
    return (log(iBudget - iWealthNext) + β * EV)

end 