##          Inner functions

#           Content:
#           1.      Next period's wealth optimsation 
#           2.      Interpolate value function
#           3.      Interpolation with the basis matrix method
#           4.      Iteration method function


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

##################################################
##          3. Basis matrix method              ##
##################################################

function        fnBasisMethodInterpolation(mOld,mNew)

    ## 1.       Numbers of elements & preallocation indices
    nOld        = length(mOld)
    nNew        = length(mNew)
    I           = zeros(Int, nNew * 2)
    J           = zeros(Int, nNew * 2)
    V           = zeros(Float64, nNew * 2)

    ## 2.       Loop 
    Count       = 0
    for (i,val) in enumerate(mNew)

        # A.    Find the old grid's index
        idx     = searchsortedlast(mOld,val)
        idx     = clamp(idx, 1, nOld - 1)

        # B.    Weights
        w       = (val - mOld[idx]) / (mOld[idx+1]-mOld[idx])

        # C.    Store weights 
        Count   += 1;   I[Count] = i;   J[Count] = idx;     V[Count] = 1 - w 
        Count   += 1;   I[Count] = i;   J[Count] = idx+1;   V[Count] = w
    end 

    ## 3.       Store results
    return  sparse(I,J,V,nNew,nOld)
end

##################################################
##          4. Iteration method function        ##
##################################################

function        fnIterationMethod(iCurrentDistribution,mPolicyWealthNext2,g)

    # A.                Unpacking
    @unpack vGridA2, vGridZ,mTransitionZ = g

    # B.                Start looping
    Res                 = zeros(size(mPolicyWealthNext2))
    for iz in eachindex(vGridZ)
        for ia in eachindex(vGridA2)

            # C.        Find weights 
            ap          = mPolicyWealthNext2[ia,iz]
            ib          = searchsortedlast(vGridA2, ap)
            ib          = clamp(ib, 1, length(vGridA2)-1)
            iu          = ib + 1
            w           = (vGridA2[iu] - ap) / (vGridA2[iu] - vGridA2[ib])
            w           = clamp(w, 0.0,1.0)

            # D.        Iterate over future labour
            Mass        = iCurrentDistribution[ia,iz]
             for izp in eachindex(vGridZ)
                Res[iu,izp] += Mass * mTransitionZ[iz,izp] * (1 - w)
                Res[ib,izp] += Mass *  mTransitionZ[iz,izp] * w
             end
        end
    end
    # E.                Done :)
    return Res
end