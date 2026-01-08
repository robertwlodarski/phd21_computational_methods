# Main functions

# 1. Aiyagari (1994) general equilibrium function

##################################################
##          1. Aiyagari (1994) GE               ##
##################################################

function fnSolveAiyagari1994(p,g)

        ## 1.           Prepare the setting
        #               Unpack parameters
        @unpack α, A, δ = p

        #               Function settings 
        iWeightOld              = 0.95
        iErrorGE                = 10.0
        iTolGE                  = 1e-8
        iTolVFI                 = 1e-8
        iTolDist                = 1e-8
        iAccInterval            = 20
        iAccStart               = 30
        iIterNumGE              = 1
        t0                      = time()
        pVerbose                = true

        #               Unpack grids 
        @unpack vGridZ, mTransitionZ, vGridA1, vGridA2 = g

        #               Preallocation 
        mValue                  = repeat(vGridA1 * 0.1,1,length(vGridZ))
        mValueNew               = zero(mValue)
        mPolicyCons             = zero(mValue)
        mPolicyWealthNext       = zero(mValue)
        iCurrentDistribution    = fill(1.0 / (length(vGridA2)*length(vGridZ)),(length(vGridA2),length(vGridZ)))
        mInterpolationConstant  = fnBasisMethodInterpolation(vGridA1,vGridA2)

        #               Prepare the stationary skill distribution 
        iErrorDistZ             = 10.0
        iTolDistZ               = 1e-8
        vDistZ                  = fill(1.0 / length(vGridZ),length(vGridZ))
        while true 
                vDistZNext      = mTransitionZ' * vDistZ
                if maximum(abs.(vDistZNext - vDistZ)) < iTolDistZ
                        break 
                end 
                vDistZ          = vDistZNext
        end
        vLabSupply              = dot(vGridZ,vDistZ)
        vMinWealth              = vGridA1[1]

        ## 2.           GE & VFI loops
        iK                      = 7.0

        #               START THE GE LOOP 
        while iErrorGE > iTolGE

                # A.    Derive the K-related terms
                iInterest       = α * A  * (iK / vLabSupply)^(α-1)-δ
                iWage           = (1 - α) * A * (iK / vLabSupply)^(α)

                # B.    Prepare to start the VFI loop 
                iErrorVFI       = 10.0
                iIterNumVFI     = 1

                #       START THE VFI LOOP
                while iErrorVFI > iTolVFI
                        #               Set the expected value
                        mExpectedValue          = mValue * mTransitionZ'
                        
                        # C.            Start iterating over a and z
                        for iz in eachindex(vGridZ)

                        # D.            Set z-related values and expected future values 
                        iLabour         = vGridZ[iz]
                        iExpValueZ      = @view mExpectedValue[:,iz]

                        for ia in eachindex(vGridA1)

                                ## 3.                   NON-ACCELERATED VERSION
                                if      (iIterNumVFI - 1) % iAccInterval == 0 || iIterNumVFI <= iAccStart    
                                        # A.            Current values
                                        a               = vGridA1[ia]
                                        iBudget         = iWage * iLabour + (1 + iInterest) * a 
                                        
                                        # B.            Optimisation
                                        iWealthNext     = fnWealthNext(iBudget,vMinWealth,iExpValueZ,p,g)
                                        iWealthNext     = max(iWealthNext, vMinWealth)
                                        iConsumption    = iBudget - iWealthNext
                                        iValue          = fnInterpolateValueFunction(iWealthNext,iBudget,iExpValueZ,p,g)

                                        # C.            Updating business
                                        mValueNew[ia,iz]        = iValue
                                        mPolicyCons[ia,iz]      = iConsumption
                                        mPolicyWealthNext[ia,iz]= iWealthNext
                                else 
                                        #               ACCELERATED VERSION
                                        # D.            Update based on the existing values
                                        a               = vGridA1[ia]
                                        iBudget         = iWage * iLabour + (1 + iInterest) * a 
                                        iConsumption    = mPolicyCons[ia,iz]
                                        iWealthNext     = mPolicyWealthNext[ia,iz]

                                        # E.            Value function
                                        iValue          = fnInterpolateValueFunction(iWealthNext,iBudget,iExpValueZ,p,g)
                                        mValueNew[ia,iz]= iValue
                                end                     

                                
                        end 
                        end 

                        # 4.            Update VFI loop items
                        iIterNumVFI     = iIterNumVFI + 1
                        iErrorVFI       = maximum(abs.(mValue - mValueNew))
                        mValue          .= mValueNew
                #       CLOSE THE VFI LOOP
                end 

                ## 5.           Wealth policy: Create a more granular grid 
                mPolicyWealthNext2      = zeros(size(iCurrentDistribution))
                if length(vGridA1) == length(vGridA2)
                        mPolicyWealthNext2 .= mPolicyWealthNext
                else 
                        mPolicyWealthNext2 = mInterpolationConstant * mPolicyWealthNext
                end 

                ## 6.           Iteration method
                iErrorDist      = 10.0
                while iErrorDist > iTolDist
                        #                       A. Compute next distribution
                        iNextDistribution       = fnIterationMethod(iCurrentDistribution,mPolicyWealthNext2,g)

                        #                       B. Compute error 
                        iErrorDist              = maximum(abs.(iNextDistribution-iCurrentDistribution))
                        iCurrentDistribution    .= iNextDistribution
                        # CLOSE THE INTERPOLATION LOOP
                end 

                ## 7.           Aggregation & updating business
                iMarginalDist   = vec(sum(iCurrentDistribution,dims=2))
                iEndoK          = dot(vGridA2,iMarginalDist)
                iErrorGE        = abs(iEndoK - iK)
                iK              = iWeightOld * iK + (1 - iWeightOld) * iEndoK

                # 7. Print cute messages 
                
                if (pVerbose && (iIterNumGE - 1) % 50 == 0) || iErrorGE < iTolGE
                        println("SUMMARY")
                        @printf("Iteration number:  %.0f \n", iIterNumGE)
                        @printf("Maximum error:     %.6f \n", iErrorGE)
                        @printf("Interest rate:     %.3f \n", iInterest)
                        @printf("Wage:              %.3f \n", iWage)
                        @printf("Elapsed time:      %.3f \n", time() - t0) 
                        flush(stdout)
                end
                iIterNumGE += 1
                #       CLOSE THE GE LOOP
        end 
return iK
end
