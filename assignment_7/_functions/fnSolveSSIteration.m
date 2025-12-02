function [Results]       = fnSolveSSIteration(Parameters,Grids)
tic; 

%% 1. Unpacking
% Parameters
pAlpha                  = Parameters.pAlpha;
pDelta                  = Parameters.pDelta;
pA                      = Parameters.pA;
pBeta                   = Parameters.pBeta;
pMu                     = Parameters.pMu;
pChi                    = Parameters.pChi;
pEta                    = Parameters.pEta;
pVerbose                = true;

% Grids
vGridA1                 = Grids.vGridA1;
vGridA2                 = Grids.vGridA2;
vGridZ                  = Grids.vGridZ;
mTransitionZ            = Grids.mTransitionZ;


%% 2. Set-up
% Iterations business: Convergence, acceleration, and GE
iWeightOld              = 0.9;
iErrorGE                = 10;
iTolGE                  = 1e-8;
iTolVFI                 = 1e-8;
iTolDist                = 1e-8;
iAccelerationInterv     = 20;
iAccelerationStart      = 30;
iIterNumGE              = 1;

% Matrices: Value function, policy functions
mValue                  = repmat(0.1*vGridA1,1,size(vGridZ,1));
mValueNew               = zeros(size(mValue));
mPolicyCons             = mValueNew;
mPolicyWealthNext       = mPolicyCons;
mPolicyLabour           = mPolicyWealthNext;
% Distributions
iCurrentDistribution    = ones(size(vGridA2,1),size(vGridZ,1)) / sum(sum(ones(size(vGridA2,1),size(vGridZ,1))));

%% 3. Stationary dist. of labour allocation
iErrorZDist             = 10;
iTolZDist               = 1e-8;
vZDist                  = ones(size(vGridZ,1),1) / sum( ones(size(vGridZ,1),1));
while iErrorZDist>iTolZDist
    vZDistNext          = mTransitionZ' * vZDist;
    iErrorZDist         = abs(vZDistNext-vZDist);
    vZDist              = vZDistNext; 
end

%% 4. GE & VFI loops

% Initial K guess
iK                      = 7.1495;
iL                      = 0.33*vGridZ'*vZDist;
iN                      = 0.33;
% Caution: iN = hours, iL = effective hours

% START GE LOOP
while iErrorGE>iTolGE
    % Derive K-related items
    iInterest           = pAlpha * pA * (iK / iL)^(pAlpha - 1) - pDelta;
    iWage               = (1 - pAlpha) * pA * (iK / iL)^(pAlpha);
    % Prepare for starting VFI loop
    iErrorVFI           = 10;
    iNumIterVFI         = 1;

    %% 4.1. START VFI LOOP
    while iErrorVFI>iTolVFI

        % Start iterating over z and a
        for zzz = 1:1:size(vGridZ,1)
            % Set z-related items
            iLabourProd                 = vGridZ(zzz);

            % Expected future values
            mExpectedValue              = mValue * mTransitionZ';
            iExpValueZ                  = mExpectedValue(:,zzz);
            iMinWealth                  = vGridA1(1);

            for aaa = 1:1:size(vGridA1,1) 

            % NON-ACCELERATED VERSION
            if (floor((iNumIterVFI-1) / iAccelerationInterv) == ((iNumIterVFI-1) / iAccelerationInterv) || iNumIterVFI <= iAccelerationStart)
            % Current values
            iWealth                     = vGridA1(aaa);
            %iBudget                     = iWage * iLabourProd + (1+iInterest) * iWealth;

            % Optimisation station
            [ap,c,n]                    = fnFastOptimisation(iWealth,iWage,iInterest,iLabourProd,iMinWealth,iExpValueZ,Parameters,Grids);
            iWealthNext                 = ap;
            iWealthNext(iWealthNext<vGridA1(1)) ...
                                        =vGridA1(1);
            iConsumption                = c;
            iLabour                     = n;

            % Interpolate value function
            iWLow                       = min(max(sum(vGridA1<iWealthNext),1),size(vGridA1,1)-1);
            iWHigh                      = iWLow+1;
            iWeightLow                  = (vGridA1(iWHigh)-iWealthNext) / (vGridA1(iWHigh)-vGridA1(iWLow));
            iExpValue                   = iWeightLow*iExpValueZ(iWLow)+(1-iWeightLow)*iExpValueZ(iWHigh);
            iValue                      = log(iConsumption)- pEta * iLabour^(1+1/pChi)/(1+1/pChi)+pBeta*iExpValue;

            % Updating
            mValueNew(aaa,zzz)          = iValue;
            mPolicyCons(aaa,zzz)        = iConsumption;
            mPolicyWealthNext(aaa,zzz)  = iWealthNext;
            mPolicyLabour(aaa,zzz)      = iLabour;

            % ACCELERATED VERSION
            else
            % Update based on the existing value
            iConsumption                = mPolicyCons(aaa,zzz);
            iWealthNext                 = mPolicyWealthNext(aaa,zzz);
            iLabour                     = mPolicyLabour(aaa,zzz);
            
            % Interpolate the value function
            iWLow                       = min(max(sum(vGridA1<iWealthNext),1),size(vGridA1,1)-1);
            iWHigh                      = iWLow+1;
            iWeightLow                  = (vGridA1(iWHigh)-iWealthNext) / (vGridA1(iWHigh)-vGridA1(iWLow));
            iExpValue                   = iWeightLow*iExpValueZ(iWLow)+(1-iWeightLow)*iExpValueZ(iWHigh);
            iValue                      = log(iConsumption)- pEta * iLabour^(1+1/pChi)/(1+1/pChi)+pBeta*iExpValue;

            % Update wealth
            mValueNew(aaa,zzz)          = iValue;
            end 
            end
        end
        % Update VFI loop items
        iNumIterVFI                 = iNumIterVFI + 1;
        iErrorVFI                   = max(abs(mValue-mValueNew),[],"all");
        mValue                      = mValueNew;  
        % END VFI LOOP
    end
    %% 4.2. Wealth policy: More granular grid
    mPolicyWealthNext2              = zeros(size(iCurrentDistribution));

    % If both are the same
    if (size(vGridA1,1)==size(vGridA2,1))
        mPolicyWealthNext2          = mPolicyWealthNext;
    else 
    end

    % If both of them aren't the same
    for zzz = 1:1:size(vGridZ,1)
        mPolicyWealthNext2(:,zzz)   = interp1(vGridA1,mPolicyWealthNext(:,zzz),vGridA2,"linear","extrap");
    end
 
    %% 4.3. Iteration method
    
    % Initialise
    iErrorDist                     = 10;
    while iErrorDist > iTolDist
    iNextDistribution               = zeros(size(iCurrentDistribution));
    % Start iterating over z and a
    for zzz = 1:1:size(vGridZ,1)
        for aaa = 1:1:size(vGridA2,1)
            % Set up
            iWealth                         = vGridA2(aaa);
            iWealthNext                     = mPolicyWealthNext2(aaa,zzz);
            iLB                             = sum(vGridA2<iWealthNext);
            iLB(iLB<=0)                     = 1;
            iLB(iLB>=size(vGridA2,1))       = size(vGridA2,1)-1;
            iUB                             = iLB + 1;
            iWeightLB                       = (vGridA2(iUB)-iWealthNext) / (vGridA2(iUB)-vGridA2(iLB));
            iWeightLB(iWeightLB<0)          = 0;
            iWeightLB(iWeightLB>1)          = 1;
            iWeightUB                       = 1 - iWeightLB;
            iMass                           = iCurrentDistribution(aaa,zzz);

            % Iterate over future labour
            for zzp = 1:1:size(vGridZ,1)
                % Upper
                iNextDistribution(iUB,zzp)  = iNextDistribution(iUB,zzp) + iMass * mTransitionZ(zzz,zzp) * iWeightUB;
                % Lower
                iNextDistribution(iLB,zzp)  = iNextDistribution(iLB,zzp) + iMass * mTransitionZ(zzz,zzp) * iWeightLB;
            end
        end
    end
    % Save
    iErrorDist          = max(abs(iNextDistribution-iCurrentDistribution),[],"all");
    iCurrentDistribution=iNextDistribution;
    end

    % Update other elements
    mPolicyLabour2      = interp1(vGridA1, mPolicyLabour, vGridA2, 'linear', 'extrap');
    iEndoL              = vGridZ' * sum(iCurrentDistribution.*mPolicyLabour2,1)';
    iEndoN              = sum(iCurrentDistribution.*mPolicyLabour2,"all");
    iMarginalDist       = sum(iCurrentDistribution,2);
    iEndoK              = vGridA2' * iMarginalDist;
    iErrorGE            = max(abs(iEndoK - iK),abs(iEndoL - iL));
    iK                  = iK .* iWeightOld + iEndoK .* (1-iWeightOld);
    iN                  = iN .* iWeightOld + iEndoN .* (1-iWeightOld);
    iL                  = iL .* iWeightOld + iEndoL .* (1-iWeightOld);
    
    
    % Print cute messages 
    if (pVerbose == true && floor((iIterNumGE-1) / 50) == ((iIterNumGE-1) / 50) || iErrorGE < iTolGE)
        fprintf('SUMMARY \n')
        fprintf('Maximum error:     %.3f \n', iErrorGE);
        fprintf('Iterest rate:      %.3f \n', iInterest);
        fprintf('Wage:              %.3f \n', iWage);
        toc;
    else 
    end
    % END GE LOOP
    iIterNumGE          = iIterNumGE + 1;
end

%% Save results
Results.mDistribution       = iCurrentDistribution;
Results.vWage               = iWage;
Results.vInterest           = iInterest;
Results.vCapitalOpt         = iK;
Results.mPolicyWealthNext   = mPolicyWealthNext2;
Results.vLabourSupply       = iL;
Results.vHours              = iN;
Results.mValue              = mValue;

end