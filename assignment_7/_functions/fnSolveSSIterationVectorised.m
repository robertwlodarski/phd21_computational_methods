function [Results]       = fnSolveSSIterationVectorised(Parameters,Grids)
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
iWeightOld              = 0.95;
iErrorGE                = 10;
iTolGE                  = 1e-5;
iTolVFI                 = 1e-5;
iTolDist                = 1e-5;
iAccelerationInterv     = 20;
iAccelerationStart      = 30;
iIterNumGE              = 1;
iWeightVFI              = 0;

% Golden ratio settings
pTolGR                  = 1e-6;
pCGR1                   = (sqrt(5)-1)/2;
pCGR2                   = (3-sqrt(5))/2;
pIterNumTolGR           = 25;

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
iK                      = 26.6 / sqrt(pEta);
iL                      = 0.94 / sqrt(pEta) *vGridZ'*vZDist;
iN                      = 0.94 / sqrt(pEta);
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


            %% A1. Golden ratio algorithm
            % Initial points
            vApLB                   = repmat(vGridA1(1),size(vGridA1,1),1);
            vApUBSafe               = max(vGridA1 *3,vGridA1+5); %PROBLEM: MASSIVE NEGATIVE VALUES 
            vApUB                   = min(repmat(vGridA1(end),size(vGridA1,1),1),vApUBSafe);
            % Initial step
            vApDiff                 = vApUB-vApLB;
            vApMid1                 = vApLB+pCGR1*vApDiff;
            vApMid2                 = vApLB+pCGR2*vApDiff;

            % Simplify the function
            fnObjective             = @(vAp,vA) fnVectorisedValues(vAp,vA,iLabourProd,iWage,iInterest,Parameters,Grids,iExpValueZ);
            vVMid1                  = fnObjective(vApMid1,vGridA1);
            vVMid2                  = fnObjective(vApMid2,vGridA1);

            % Run the golden ratio algorithm
            for iii = 1:1:25
                % Identify better region
                iIndexRegion            = (vVMid2 >= vVMid1);
                vApLB(iIndexRegion)     = vApMid1(iIndexRegion);
                vApUB(~iIndexRegion)    = vApMid2(~iIndexRegion);

                % Update points [if vVMid1 >= vVmid2]
                vApDiff(iIndexRegion)   = vApUB(iIndexRegion)-vApLB(iIndexRegion);
                vApMid1(iIndexRegion)   = vApMid2(iIndexRegion);
                vApMid2(iIndexRegion)   = vApLB(iIndexRegion)+pCGR2*vApDiff(iIndexRegion);
                vVMid1(iIndexRegion)    = vVMid2(iIndexRegion);
                vVMid2(iIndexRegion)    = fnObjective(vApMid2(iIndexRegion),vGridA1(iIndexRegion));

                % Update points [if vVMid1 < vVmid2]
                vApDiff(~iIndexRegion)  = vApUB(~iIndexRegion)-vApLB(~iIndexRegion);
                vApMid2(~iIndexRegion)  = vApMid1(~iIndexRegion);
                vApMid1(~iIndexRegion)  = vApLB(~iIndexRegion)+pCGR1*vApDiff(~iIndexRegion);
                vVMid2(~iIndexRegion)   = vVMid1(~iIndexRegion);
                vVMid1(~iIndexRegion)   = fnObjective(vApMid1(~iIndexRegion),vGridA1(~iIndexRegion));
            end

            % Save results
            vAp                         = (vApLB + vApUB) / 2;
            [vV,vC,vN]                  = fnObjective(vAp,vGridA1);

            % Updating
            mValueNew(:,zzz)            = vV;
            mPolicyCons(:,zzz)          = vC;
            mPolicyWealthNext(:,zzz)    = vAp;
            mPolicyLabour(:,zzz)        = vN;

        end
        % Update VFI loop items
        iNumIterVFI                 = iNumIterVFI + 1;
        iErrorVFI                   = max(abs(mValue-mValueNew),[],"all");
        mValue                      = iWeightVFI*mValue + (1-iWeightVFI)*mValueNew;
        fprintf("Error:  %.10f\n",iErrorVFI);
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
        fprintf('======================================== \n')
        fprintf('             Iteration %.0f              \n',iIterNumGE);
        fprintf('======================================== \n')
        fprintf('SUMMARY \n')
        fprintf('Maximum error:     %.3f \n', iErrorGE);
        fprintf('Iterest rate:      %.3f \n', iInterest);
        fprintf('Wage:              %.3f \n', iWage);
        fprintf('IMPLIED VALUES \n')
        fprintf('K:                 %.3f \n', iEndoK);
        fprintf('L:                 %.3f \n', iEndoL);
        fprintf('N:                 %.3f \n', iEndoN);
        fprintf('NEW VALUES \n')
        fprintf('K:                 %.3f \n', iK);
        fprintf('L:                 %.3f \n', iL);
        fprintf('N:                 %.3f \n', iN);
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