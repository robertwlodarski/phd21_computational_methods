function [Results]          = fnAiyagari1994Transition(Start,End,Parameters,Grids)

%% 1. Unpacking
% Parameters
ParametersUsed          = Parameters;
pT                      = Parameters.pT;
pAlpha                  = ParametersUsed.pAlpha;
pDelta                  = ParametersUsed.pDelta;
pBeta                   = ParametersUsed.pBeta;

% Grids
vGridA1                 = Grids.vGridA1;
vGridA2                 = Grids.vGridA2;
vGridZ                  = Grids.vGridZ;
mTransitionZ            = Grids.mTransitionZ;
vAPath                  = linspace(Start,End,pT)';   
%% 2. Compute the initial and terminal equilibria

% Terminal state
ParametersUsed.pA       = End;
ResultsT                = fnSolveAiyagari1994Iteration(ParametersUsed,Grids);

% Initial state
ParametersUsed.pA       = Start;
Results0                = fnSolveAiyagari1994Iteration(ParametersUsed,Grids);

%% 3. Set-up
% Iterations business: Convergence, acceleration, and GE
iWeightOld              = 0.9;
iErrorK                 = 10;
iTolK                   = 1e-4;
iIterNumK               = 1;

%% 4. Stationary dist. of labour allocation
iErrorZDist             = 10;
iTolZDist               = 1e-8;
vZDist                  = ones(size(vGridZ,1),1) / sum( ones(size(vGridZ,1),1));
while iErrorZDist>iTolZDist
    vZDistNext          = mTransitionZ' * vZDist;
    iErrorZDist         = abs(vZDistNext-vZDist);
    vZDist              = vZDistNext; 
end
vLabSupply              = vGridZ'*vZDist;

%% 5. TRANSITIONAL DYNAMICS

% Capital guess
iK                      = linspace(Results0.vCapitalOpt,ResultsT.vCapitalOpt,pT)';
iEndoK                  = zeros(size(iK));
iGini                   = zeros(size(iK));

% Initialise values
iValueNext              = ResultsT.mValue;
mValueNew               = zeros(size(iValueNext));
mPolicyWealth           = zeros([pT,size(vGridA2,1),size(vGridZ,1)]);

while iErrorK > iTolK
    % Initialise
    iCurrentDistribution    = Results0.mDistribution;
    iValueNext              = ResultsT.mValue;

    % Compute endogenous objects
    iInterest               = pAlpha .* vAPath .* (iK / vLabSupply).^(pAlpha - 1) - pDelta;
    iWage                   = (1 - pAlpha) .* vAPath .* (iK / vLabSupply).^(pAlpha);

    %% 5.1. Value & policy functions
    for ttt = pT:(-1):1
        
        % A. Matrices: Value function, policy functions
        mPolicyWealthNext       = zeros(size(iValueNext));

        % B. Run the loop
        % Start iterating over z and a
        for zzz = 1:1:size(vGridZ,1)
            % Set z-related items
            iLabour                         = vGridZ(zzz);

            % Expected future values
            mExpectedValue                  = iValueNext * mTransitionZ';
            iExpValueZ                      = mExpectedValue(:,zzz);
            iMinWealth                      = vGridA1(1);

            for aaa = 1:1:size(vGridA1,1)

                % Current values
                iWealth                     = vGridA1(aaa);
                iBudget                     = iWage(ttt) * iLabour + (1+iInterest(ttt)) * iWealth;

                % Optimisation station
                iWealthNext                 = fnFastOptimisation(iBudget,iMinWealth,iExpValueZ,ParametersUsed,Grids);
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

        % C. Update the wealth grid
        mPolicyWealthNext2                  = zeros(size(iCurrentDistribution));

        % If both are the same
        if (size(vGridA1,1)==size(vGridA2,1))
            mPolicyWealthNext2              = mPolicyWealthNext;
        else
            % If both of them aren't the same
            for zzz = 1:1:size(vGridZ,1)
                mPolicyWealthNext2(:,zzz)   = interp1(vGridA1,mPolicyWealthNext(:,zzz),vGridA2,"linear","extrap");
            end
        end

        % D. Save the results
        mPolicyWealth(ttt,:,:)              = mPolicyWealthNext2;
        iValueNext                          = mValueNew;
    end

    %% 5.2. Distributions (iteration method) and aggregate to obtain new K
    for ttt = 1:1:pT
        
        % A. Finding distributions
        iNextDistribution                       = zeros(size(iCurrentDistribution));
        % Start iterating over z and a
        for zzz = 1:1:size(vGridZ,1)
            for aaa = 1:1:size(vGridA2,1)
                % Set up
                iWealthNext                     = mPolicyWealth(ttt,aaa,zzz);
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
        iCurrentDistribution=iNextDistribution;

        % B. Aggregation station
        iMarginalDist       = sum(iCurrentDistribution,2);
        iEndoK(ttt)         = vGridA2' * iMarginalDist;
        iGini(ttt)          = 1 / (2*iEndoK(ttt)) * sum((iMarginalDist * iMarginalDist') .* abs( repmat(vGridA2',size(vGridA2,1),1)-repmat(vGridA2,1,size(vGridA2,1))),'all');
    end 

    % C. Update
    iErrorK                 = max(abs(iEndoK - iK),[],"all");
    iK                      = iK .* iWeightOld + iEndoK .* (1-iWeightOld);

    % Print error and update iteration count
    fprintf("Iteration: %.0f, K Error: %.4f \n", iIterNumK, iErrorK);
    iIterNumK               = iIterNumK + 1;

    % End the equilibrium loop
end

Results.vK                  = iK;
Results.vGrini              = iGini;
Results.vK0                 = Results0.vCapitalOpt;
Results.vKT                 = ResultsT.vCapitalOpt;

end 