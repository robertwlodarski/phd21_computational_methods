function [EndoK, DistributionImplied, ValueFunction] =fnIntermediateEquilibriumStep(FutureValue,GuessK,FutureDistribution,vLabSupply, ParametersUsed,Grids)

% ARCHIVED, NOT TO BE USED. 

%% 1. Unpacking
% Parameters
pAlpha                  = ParametersUsed.pAlpha;
pDelta                  = ParametersUsed.pDelta;
pA                      = ParametersUsed.pA;
pBeta                   = ParametersUsed.pBeta;

% Grids
vGridA1                 = Grids.vGridA1;
vGridA2                 = Grids.vGridA2;
vGridZ                  = Grids.vGridZ;
mTransitionZ            = Grids.mTransitionZ;

%% 2. Pre-iteration elements
% Compute endogenous objects
iInterest               = pAlpha * pA * (GuessK / vLabSupply)^(pAlpha - 1) - pDelta;
iWage                   = (1 - pAlpha) * pA * (GuessK / vLabSupply)^(pAlpha);

% Matrices: Value function, policy functions
mValueNew               = zeros(size(FutureValue));
mPolicyWealthNext       = mValueNew;
% Distribution guess (starting from future for speed)
iCurrentDistribution    = FutureDistribution;

%% 3. Value function

% Start iterating over z and a
for zzz = 1:1:size(vGridZ,1)
    % Set z-related items
    iLabour                         = vGridZ(zzz);

    % Expected future values
    mExpectedValue                  = FutureValue * mTransitionZ';
    iExpValueZ                      = mExpectedValue(:,zzz);
    iMinWealth                      = vGridA1(1);

    for aaa = 1:1:size(vGridA1,1)

        % Current values
        iWealth                     = vGridA1(aaa);
        iBudget                     = iWage * iLabour + (1+iInterest) * iWealth;

        % Optimisation station
        iWealthNext                 = fnFastOptimisation(iBudget,iMinWealth,iExpValueZ,Parameters,Grids);
        iWealthNext(iWealthNext<vGridA1(1)) ...
            =vGridA1(1);
        iConsumption                = iBudget-iWealthNext;

        % Interpolate value function
        iWLow                       = min(max(sum(vGridA1<iWealthNext),1),size(vGridA1,1)-1);
        iWHigh                      = iWLow+1;
        iWeightLow                  = (vGridA1(iWHigh)-iWealthNext) / (vGridA1(iWHigh)-vGridA1(iWLow));
        iExpValue                   = iWeightLow*iExpValueZ(iWLow)+(1-iWeightLow)*iExpValueZ(iWHigh);
        iValue                      = log(iConsumption)+pBeta*iExpValue;

        % Updating
        mValueNew(aaa,zzz)          = iValue;
        mPolicyWealthNext(aaa,zzz)  = iWealthNext;
    end
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
iErrorDist                      = 10;
iTolDist                        = 1e-8;
while iErrorDist > iTolDist
    iNextDistribution               = zeros(size(iCurrentDistribution));
    % Start iterating over z and a
    for zzz = 1:1:size(vGridZ,1)
        for aaa = 1:1:size(vGridA2,1)
            % Set up
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
iMarginalDist       = sum(iCurrentDistribution,2);
iEndoK              = vGridA2' * iMarginalDist;

%% RESULTS
EndoK               = iEndoK;
DistributionImplied = iCurrentDistribution;
ValueFunction       = mValueNew;

end 