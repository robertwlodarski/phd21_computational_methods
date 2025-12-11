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

%% 2. Stationary dist. of labour allocation
iErrorZDist             = 10;
iTolZDist               = 1e-8;
vZDist                  = ones(size(vGridZ,1),1) / sum( ones(size(vGridZ,1),1));
while iErrorZDist>iTolZDist
    vZDistNext          = mTransitionZ' * vZDist;
    iErrorZDist         = abs(vZDistNext-vZDist);
    vZDist              = vZDistNext; 
end

%% 3. Set-up
% Iterations business: Convergence, acceleration, and GE
iWeightKtoL             = 0.9999;
iErrorGE                = 10;
iErrorPFI               = 10;
iTolGE                  = 1e-5;
iTolPFI                 = 1e-4;
iTolDist                = 1e-5;
iIterNumGE              = 1;
iWeightPolicy           = 0.95;

% Initial K guess
iKtoL                   = ( (pAlpha * pA) / (1 / pBeta -1 + pDelta) ) ^ (1 / (1 - pAlpha));

% Other initialised elements
mPolicyWealthNext       = repmat(vGridA1,1,size(vGridZ,1))+repmat(vGridZ',size(vGridA1,1),1);
mPolicyCons             = zeros(size(mPolicyWealthNext));
mPolicyLabour           = mPolicyWealthNext;
% Distributions
iCurrentDistribution    = ones(size(vGridA2,1),size(vGridZ,1)) / sum(sum(ones(size(vGridA2,1),size(vGridZ,1))));


%% 4. GE & PFI loops

% START GE LOOP
while (iErrorGE>iTolGE && iErrorPFI > iTolPFI)
    % Derive K-related items
    iInterest       = pAlpha * pA * (iKtoL)^(pAlpha - 1) - pDelta;
    iWage           = (1 - pAlpha) * pA * (iKtoL)^(pAlpha);

    %% 4.1. START PFI LOOP
    %while iErrorPFI > iTolPFI
       
    %% A. Interpolate mApp
    mA              = repmat(vGridA1,1,size(vGridZ,1));
    mAp             = mPolicyWealthNext;
    mApp            = zeros(size(mPolicyWealthNext));

    % Loop over the labout productivity grid
    for zzz = 1:1:size(vGridZ,1)
        mApp(:,zzz)         = interp1(vGridA1,mPolicyWealthNext(:,zzz),mPolicyWealthNext(:,zzz),'linear','extrap');
    end

    %% B. Use (Budget) and (Labour) to back out c'
    % Remember that you're one period ahead
    mPsip           = pMu / 2 * (mApp.^2 ./ mAp - 2* mApp + mAp);
    mAAA            = (1+iInterest)*mAp - mApp - mPsip;
    mZ              = repmat(vGridZ',size(vGridA1,1),1);
    mBBB            = (iWage * mZ).^2 / pEta;
    mCp             = (sqrt(mAAA.^2 + 4* mBBB)+mAAA)/2;

    %% C. Get current optimal c
    mPsi1           = pMu * (mAp ./ mA -1);
    mDenominator    = 1 + mPsi1;
    mPsi2p          = pMu / 2 * (- (mApp./mAp).^2 + 1);
    mInnerTerm      = mCp.^(-1).*(1+iInterest-mPsi2p);
    mNumerator      = pBeta * mInnerTerm * mTransitionZ';
    mC              = (mNumerator ./ mDenominator).^(-1);

    %% D. New labour (things can only get better?)
    mN              = iWage * mZ ./ (mC * pEta);

    %% E. New policy
    mPsi            = pMu / 2 * (mAp.^2 ./ mA - 2* mAp + mA);
    mBudget         = (1+iInterest)*mA + iWage * mZ .* mN;
    mApNew          = max(mBudget-mPsi-mC,vGridA1(1));

    %% F. Update consumption and labour
    for zzz = 1:1:size(vGridZ,1)
        mApp(:,zzz)         = interp1(vGridA1,mApNew(:,zzz),mApNew(:,zzz),'linear','extrap');
    end
    mPsi1           = pMu * (mApNew ./ mA -1);
    mDenominator    = 1 + mPsi1;
    mPsi2p          = pMu / 2 * (- (mApp./mApNew).^2 + 1);
    mInnerTerm      = mCp.^(-1).*(1+iInterest-mPsi2p);
    mNumerator      = pBeta * mInnerTerm * mTransitionZ';
    mC              = (mNumerator ./ mDenominator).^(-1);
    mN              = iWage * mZ ./ (mC * pEta);
    

    %fprintf("PFI error %.8f \n",iErrorPFI);
    % Finish PFI here
    %end

    %% 4.2. Wealth policy: More granular grid
    mPolicyWealthNext2              = zeros(size(iCurrentDistribution));

    % If both are the same
    if (size(vGridA1,1)==size(vGridA2,1))
        mPolicyWealthNext2          = mApNew;
    else 
    end

    % If both of them aren't the same
    for zzz = 1:1:size(vGridZ,1)
        mPolicyWealthNext2(:,zzz)   = interp1(vGridA1,mApNew(:,zzz),vGridA2,"linear","extrap");
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

    %% 4.4. Update & compute error
    iErrorPFI           = max(abs(mApNew - mAp),[],"all");
    mPolicyLabour       = mN;


    % Update other elements
    mPolicyLabour2      = interp1(vGridA1, mPolicyLabour, vGridA2, 'linear', 'extrap');
    iEndoL              = vGridZ' * sum(iCurrentDistribution.*mPolicyLabour2,1)';
    iEndoN              = sum(iCurrentDistribution.*mPolicyLabour2,"all");
    iMarginalDist       = sum(iCurrentDistribution,2);
    iEndoK              = vGridA2' * iMarginalDist;
    iEndoKtoL           = iEndoK/iEndoL;
    iDiff               = iEndoKtoL- iKtoL;
    iErrorGE            = max(abs(iDiff));

    % Update K and policies
    mPolicyWealthNext   = iWeightPolicy .* mAp + (1-iWeightPolicy).*mApNew; 
    iKtoL               = iWeightKtoL * iKtoL + (1-iWeightKtoL)*iEndoKtoL;

    
    % Print cute messages 
    if (pVerbose == true && floor((iIterNumGE-1) / 100) == ((iIterNumGE-1) / 100) || iErrorGE < iTolGE)
        fprintf('======================================== \n')
        fprintf('             Iteration %.0f              \n',iIterNumGE);
        fprintf('======================================== \n')
        fprintf('SUMMARY \n')
        fprintf('GE error:          %.3f \n', iErrorGE);
        fprintf('PFI error:         %.3f \n', iErrorPFI);
        fprintf('Iterest rate:      %.3f \n', iInterest);
        fprintf('Wage:              %.3f \n', iWage);
        fprintf('IMPLIED VALUES \n')
        fprintf('K-to-L:            %.3f \n', iEndoKtoL);
        fprintf('NEW VALUES \n')
        fprintf('K-to-L:            %.3f \n', iKtoL);
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
Results.vCapitalOpt         = iEndoK;
Results.mPolicyWealthNext   = mPolicyWealthNext2;
Results.vLabourSupply       = iEndoL;
Results.vHours              = iEndoN;
Results.mC                  = mC;
Results.mN                  = mN;
%Results.mValue              = mValue;

save("_results/DistributionSS","iCurrentDistribution");
save("_results/Policy1SS","mApNew");
save("_results/Policy2SS","mPolicyWealthNext2");
save("_results/KtoLSS","iKtoL");

end


% ARCHIVE: VFI LOOP WTH GOLDEN SEARCH AND HOWARD ACCELERATOR [NOT WORKING]
% % Golden ratio settings
% pTolGR                  = 1e-8;
% pCGR1                   = (sqrt(5)-1)/2;
% pCGR2                   = (3-sqrt(5))/2;
% pIterNumTolGR           = 25;
% while iErrorVFI>iTolVFI
%
%     % Start iterating over z and a
%     for zzz = 1:1:size(vGridZ,1)
%         % Set z-related items
%         iLabourProd                 = vGridZ(zzz);
%
%         % Expected future values
%         mExpectedValue              = mValue * mTransitionZ';
%         iExpValueZ                  = mExpectedValue(:,zzz);
%         iMinWealth                  = vGridA1(1);
%
%         %% A. POLICY FUNCTIONS
%         % A1. Golden ratio algorithm
%         % Initial points
%         vApLB                   = repmat(vGridA1(1),size(vGridA1,1),1);
%         vApUBSafe               = (1+iInterest)*vGridA1 + iWage * iLabourProd * 5;
%         vApUB                   = min(repmat(vGridA1(end),size(vGridA1,1),1),vApUBSafe);
%
%         % Prohibit radical jumps for the Golden ratio algorithm
%         if iNumIterVFI >1
%             % For lower bound
%             vStepL              = 2;
%             vApLB               = max(vApLB,mPolicyWealthNext(:,zzz)-vStepL);
%             % For upper bound
%             vStepU              = 5;
%             vApUB               = min(vApUB,mPolicyWealthNext(:,zzz)+vStepU);
%             % Safety check
%             vApUB               = max(vApLB+1e-5,vApUB);
%         end
%
%         % Initial step
%         vApDiff                 = vApUB-vApLB;
%         vApMid1                 = vApLB+pCGR1*vApDiff;
%         vApMid2                 = vApLB+pCGR2*vApDiff;
%
%         % Simplify the function
%         fnObjective             = @(vAp,vA) fnVectorisedValues(vAp,vA,iLabourProd,iWage,iInterest,Parameters,Grids,iExpValueZ);
%         vVMid1                  = fnObjective(vApMid1,vGridA1);
%         vVMid2                  = fnObjective(vApMid2,vGridA1);
%
%         % Run the golden ratio algorithm
%         for iii = 1:1:pIterNumTolGR
%             % Identify better region
%             iIndexRegion            = (vVMid2 >= vVMid1);
%             vApLB(iIndexRegion)     = vApMid1(iIndexRegion);
%             vApUB(~iIndexRegion)    = vApMid2(~iIndexRegion);
%
%             % Update points [if vVMid1 >= vVmid2]
%             vApDiff(iIndexRegion)   = vApUB(iIndexRegion)-vApLB(iIndexRegion);
%             vApMid1(iIndexRegion)   = vApMid2(iIndexRegion);
%             vApMid2(iIndexRegion)   = vApLB(iIndexRegion)+pCGR2*vApDiff(iIndexRegion);
%             vVMid1(iIndexRegion)    = vVMid2(iIndexRegion);
%             vVMid2(iIndexRegion)    = fnObjective(vApMid2(iIndexRegion),vGridA1(iIndexRegion));
%
%             % Update points [if vVMid1 < vVmid2]
%             vApDiff(~iIndexRegion)  = vApUB(~iIndexRegion)-vApLB(~iIndexRegion);
%             vApMid2(~iIndexRegion)  = vApMid1(~iIndexRegion);
%             vApMid1(~iIndexRegion)  = vApLB(~iIndexRegion)+pCGR1*vApDiff(~iIndexRegion);
%             vVMid2(~iIndexRegion)   = vVMid1(~iIndexRegion);
%             vVMid1(~iIndexRegion)   = fnObjective(vApMid1(~iIndexRegion),vGridA1(~iIndexRegion));
%         end
%
%         % Save results from the optimisation step
%         vAp                         = (vApLB + vApUB) / 2;
%         [vV,vC,vN]                  = fnObjective(vAp,vGridA1);
%
%         % Updating
%         mValueNew(:,zzz)            = vV;
%         mPolicyCons(:,zzz)          = vC;
%         mPolicyWealthNext(:,zzz)    = vAp;
%         mPolicyLabour(:,zzz)        = vN;
%
%     end
%
%     %% B. HOWARD ACCELERATOR
%
%     % Compute the flow utility component
%     mFlowUtility                    = log(mPolicyCons) - pEta * mPolicyLabour .^(1 + 1 / pChi) / (1 + 1 / pChi);
%
%     % Run the accelerated version
%     if iNumIterVFI > iAccelerationStart % Let it solve properly a few times
%         for iii = 1:1:iAccelerationInterv
%             % Update the expected value
%             mExpectedValue                      = mValueNew * mTransitionZ';
%
%             % Compute the ever-changing continuation value
%             mContinuationValue                  = zeros(size(mValue));
%             for zzz = 1:1:size(vGridZ,1)
%                 mContinuationValue(:,zzz)       = interp1(vGridA1,mExpectedValue(:,zzz),mPolicyWealthNext(:,zzz),'linear','extrap');
%             end
%             mValueNew                           = mFlowUtility + pBeta * mContinuationValue;
%         end
%     end
%
%     %% C. UPDATING BUSINESS
%
%     % Update VFI loop items
%     iNumIterVFI                 = iNumIterVFI + 1;
%     iErrorVFI                   = max(abs(mValue-mValueNew),[],"all");
%     mValue                      = iWeightVFI*mValue + (1-iWeightVFI)*mValueNew;
%
%     if (pVerbose == true && floor((iNumIterVFI-1) / 50) == ((iNumIterVFI-1) / 50) || iErrorVFI < iTolVFI)
%     fprintf('--------------------------- \n')
%     fprintf('Iteration:         %.0f \n', iNumIterVFI);
%     fprintf('VFI error:         %.3f \n', iErrorVFI);
%     else
%     end
% end

% ARCHIVE OLD VALUE COMPUTATIONS
% Save the next value function
% mAp             = mPolicyWealthNext2;
% mA              = repmat(vGridA2,1,size(vGridZ,1));
% mPsi            = pMu / 2 * (mAp.^2 ./ mA - 2* mAp + mA);
% mAAA            = (1+iInterest)*mA - mAp - mPsi;
% mZ              = repmat(vGridZ',size(vGridA2,1),1);
% mBBB            = (iWage * mZ).^2 / pEta;
% mC              = (sqrt(mAAA.^2 + 4* mBBB)+mAAA)/2;
% mN              = iWage * mZ ./ (mC * pEta);
% mValue          = log(mC) - pEta * mN.^(1+1/pChi)/(1+1/pChi)


% The slope method
% % Slope method
% if iIterNumGE==1
%     iKtoLNext           = iKtoL .* iWeightOld + iEndoKtoL .* (1-iWeightOld);
% else
%     % Computationss
%     dError              = iDiff - iDiffOld;
%     dStep               = iKtoL - iKtoLOld;
%     if abs(dStep)<1e-5
%         iSlope          = 1;
%     else
%         iSlope          = (dError) / (dStep);
%     end
%     iUpdate             = -iDiff / iSlope;
%     iMaxChange          = 0.2*iKtoL;
%     iUpdate             = max(-iMaxChange,min(iMaxChange,iUpdate));
%     iKtoLNext           = iKtoL + iUpdate;
% end
% iKtoLOld                = iKtoL;
% iDiffOld                = iDiff;
% iKtoL                   = max(0.01,iKtoLNext);