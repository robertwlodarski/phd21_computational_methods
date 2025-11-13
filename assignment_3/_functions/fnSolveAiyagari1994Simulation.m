function [Results]       = fnSolveAiyagari1994Simulation(Parameters,Grids)
tic; 

%% 1. Unpacking
% Parameters
pAlpha                  = Parameters.pAlpha;
pDelta                  = Parameters.pDelta;
pA                      = Parameters.pA;
pBeta                   = Parameters.pBeta;
pVerbose                = true;

% Grids
vGridA1                 = Grids.vGridA1;
vGridA2                 = Grids.vGridA2;
vGridZ                  = Grids.vGridZ;
mTransitionZ            = Grids.mTransitionZ;

% Simulation items
pT                      = Parameters.pT;
pN                      = Parameters.pN;
rng(1997,"twister");
mLabourShocks           = rand(pT,pN);
mLabour                 = zeros(size(mLabourShocks));
mAssets                 = mLabour;
mAssetsNext             = mAssets;
mConsumption            = mAssetsNext;
mEarnings               = mConsumption;
% Assumption: Everyone starts w/ tiny assets.
mAssets(1,:)            = vGridA1(2);

% Generate labour shocks
mLabour                 = fnGenerateProductivityValues(mLabourShocks, Grids);

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

%% 3. Stationary dist. of labour allocation
iErrorZDist             = 10;
iTolZDist               = 1e-8;
vZDist                  = ones(size(vGridZ,1),1) / sum( ones(size(vGridZ,1),1));
while iErrorZDist>iTolZDist
    vZDistNext          = mTransitionZ' * vZDist;
    iErrorZDist         = abs(vZDistNext-vZDist);
    vZDist              = vZDistNext; 
end
vLabSupply              = vGridZ'*vZDist;

%% 4. VFI & GE


% Initial K guess
iK                      = 7.1495;

% START GE LOOP
while iErrorGE>iTolGE

    % Derive K-related items
    iInterest           = pAlpha * pA * (iK / vLabSupply)^(pAlpha - 1) - pDelta;
    iWage               = (1 - pAlpha) * pA * (iK / vLabSupply)^(pAlpha);
    % Prepare for starting VFI loop
    iErrorVFI           = 10;
    iNumIterVFI         = 1;


    % Prepare for starting VFI loop
    iErrorVFI           = 10;
    iNumIterVFI         = 1;

    %% 4.1. START VFI LOOP
    while iErrorVFI>iTolVFI

        % Start iterating over z and a
        for zzz = 1:1:size(vGridZ,1)
            % Set z-related items
            iLabour                     = vGridZ(zzz);

            % Expected future values
            mExpectedValue              = mValue * mTransitionZ';
            iExpValueZ                  = mExpectedValue(:,zzz);
            iMinWealth                  = vGridA1(1);

            for aaa = 1:1:size(vGridA1,1)

                % NON-ACCELERATED VERSION
                if (floor((iNumIterVFI-1) / iAccelerationInterv) == ((iNumIterVFI-1) / iAccelerationInterv) || iNumIterVFI <= iAccelerationStart)
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
                    mPolicyCons(aaa,zzz)        = iConsumption;
                    mPolicyWealthNext(aaa,zzz)  = iWealthNext;

                    % ACCELERATED VERSION
                else
                    % Update based on the existing value
                    iConsumption                = mPolicyCons(aaa,zzz);
                    iWealthNext                 = mPolicyWealthNext(aaa,zzz);

                    % Interpolate the value function
                    iWLow                       = min(max(sum(vGridA1<iWealthNext),1),size(vGridA1,1)-1);
                    iWHigh                      = iWLow+1;
                    iWeightLow                  = (vGridA1(iWHigh)-iWealthNext) / (vGridA1(iWHigh)-vGridA1(iWLow));
                    iExpValue                   = iWeightLow*iExpValueZ(iWLow)+(1-iWeightLow)*iExpValueZ(iWHigh);
                    iValue                      = log(iConsumption)+pBeta*iExpValue;

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
    
    
    for ttt= 1:1:pT
        % Preload assets from the previous period
        if ttt > 1
            mAssets(ttt,:)  = mAssetsNext(ttt-1,:);
        else
        end
    
        % Derive the current distribution
        vLabour             = mLabour(ttt,:)';
        vAssets             = mAssets(ttt,:)';
        iCurrentDistribution= fnDeriveStochasticDistribution(vAssets,vLabour,pN,Grids);

        %% 4.2. SIMULATIONS

        % Start loop for each agents
        for nnn = 1:1:pN

            % Load & compute values
            z_index                     = mLabour(ttt,nnn);
            z                           = vGridZ(z_index);
            a                           = mAssets(ttt,nnn);
            BudgetConstraint            = iWage * z + a * (1 + iInterest);
            vAssetsNextOptions          = mPolicyWealthNext(:,z_index);
            a_next                      = interp1(vGridA1,vAssetsNextOptions,a,"linear");

            % Save agent's results
            mAssetsNext(ttt,nnn)        = a_next;
            mConsumption(ttt,nnn)       = BudgetConstraint-a_next;
            mEarnings(ttt,nnn)          = z*iWage;
        end

        % Update the distribution
    end

    % Update other elements
    iMarginalDist       = sum(iCurrentDistribution,2);
    iEndoK              = vGridA2' * iMarginalDist;
    iErrorGE            = abs(iEndoK - iK);
    iK                  = iK .* iWeightOld + iEndoK .* (1-iWeightOld);

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
Results.mAssets                     = mAssets;
Results.mConsumption                = mConsumption;
Results.mAssetsNext                 = mAssetsNext;
Results.vWage                       = iWage;
Results.vInterest                   = iInterest;
Results.vCapitalOpt                 = iK;
Results.mPolicyWealthNext           = mPolicyWealthNext;
Results.vLabourSupply               = vLabSupply;

end