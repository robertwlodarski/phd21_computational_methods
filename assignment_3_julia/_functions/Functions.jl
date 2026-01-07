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

        #               Unpack grids 
        @unpack vGridZ, mTransitionZ, vGridA1, vGridA2 = g

        #               Preallocation 
        mValue                  = repeat(vGridA1 * 0.1,1,length(vGridZ))
        mValueNew               = zero(mValue)
        mPolicyCons             = zero(mValue)
        mPolicyWealthNext       = zero(mValue)
        iCurrentDistribution    = fill(1.0 / (length(vGridA2)*length(vGridZ)),(length(vGridA2),length(vGridZ)))
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
                mExpectedValue          = mValue * mTransitionZ'

                # B.    Prepare to start the VFI loop 
                iErrorVFI       = 10.0
                iIterNumVFI     = 1

                #       START THE VFI LOOP
                while iErrorVFI > iTolVFI
                        
                        # C.            Start iterating over a and z
                        for iz in eachindex(vGridZ)

                        # D.            Set z-related values and expected future values 
                        iLabour         = vGridZ[iz]
                        iExpValueZ      = @view mExpectedValue[:,iz]

                        for ia in eachindex(vGridA1)

                                ## 3.   NON-ACCELERATED VERSION
                                if      (iIterNumVFI - 1) % iAccInterval == 0 || iIterNumVFI <= iAccStart    
                                        # A.            Current values
                                        a               = vGridA1[ia]
                                        iBudget         = iWage * iLabour + (1 + iInterest) * a 
                                        
                                        # B.            Optimisation
                                        iWealthNext     = fnWealthNext(iBudget,vMinWealth,iExpValueZ,p,g)
                                else 
                                end                     

                                #       ACCELERATED VERSION
                        end 
                        end 

                #       CLOSE THE VFI LOOP
                end 
                #       CLOSE THE GE LOOP
        end 

end
