function [Simulations]  = fnGenerateSimulations(N,T,Results, Parameters, Grids)
    %% Load in results
    mDistribution       = Results.mDistribution;
    mPolicyWealthNext   = Results.mPolicyWealthNext;
    vWage               = Results.vWage;
    vInterest           = Results.vInterest;


    %% Generate draw from the equilibrium distribution for T =  1 
    NumStates           = numel(mDistribution);
    vFlatDistribution   = mDistribution(:);
    vStateIndices       = 1:NumStates;
    vSampleIndices      = randsample(vStateIndices,N, true, vFlatDistribution);
    [a_ind,z_ind]       = ind2sub(size(mDistribution),vSampleIndices);
    
    %% Simulate agents
    % Prepare space
    mWealthIndices      = zeros(T,N);
    mLabourIndices      = mWealthIndices;
    mWealthIndices(1,:) = a_ind;
    mLabourIndices(1,:) = z_ind;

    % Grid elements
    vGridA2             = Grids.vGridA2;
    vGridZ              = Grids.vGridZ;
    pNA                 = size(vGridA2,1);
    pNZ                 = size(vGridZ,1);
    mTransitionZ        = Grids.mTransitionZ;

    % Labour shocks
    rng(1997,"twister");
    mLabourShocks       = rand(T,N);
    mCumTransitionZ     = cumsum(mTransitionZ,2);
    mBernoulliShocks    = rand(T,N);
    
    % Simulate further peridods
    for ttt=2:1:T
        % Previous wealth & labour
        iIndexPreva                 = mWealthIndices(ttt-1,:);
        iIndexPrevz                 = mLabourIndices(ttt-1,:);
        iIndices                    = sub2ind([pNA,pNZ],iIndexPreva,iIndexPrevz);

        % Generate the wealth from policy
        iNextWealth                 = mPolicyWealthNext(iIndices);
        iLB                         = sum(vGridA2 <= iNextWealth);
        iLB(iLB<=0)                 = 1;
        iLB(iLB >= size(vGridA2,1)) = size(vGridA2,1)-1;
        iUB                         = iLB + 1;
        iRatio                      = (iNextWealth-vGridA2(iLB)')./(vGridA2(iUB)'-vGridA2(iLB)');
        iBernoulli                  = iRatio>mBernoulliShocks(ttt,:);
        iIndicesNW                  = iLB + iBernoulli;

        % Generate new wealth 
        mWealthIndices(ttt,:)   = iIndicesNW;
        
        % New labour indices
        for iii=1:1:N
            mLabourIndices(ttt,iii) = sum(mLabourShocks(ttt,iii) > mCumTransitionZ(mLabourIndices(ttt-1,iii),:))+1;
        end
    end

    % Save results for wealth and labour
    Simulations.mWealth     = vGridA2(mWealthIndices);
    Simulations.mLabour     = vGridZ(mLabourIndices);
    Simulations.mWages      = Simulations.mLabour .* vWage;
    Simulations.mEarnings   = Simulations.mWages + (1+vInterest).*Simulations.mWealth;

    % Consumption (dirty, without last period)
    Simulations.mConsumption            = NaN(T,N);
    Simulations.mConsumption(1:end-1,:) = Simulations.mEarnings(1:end-1,:)-Simulations.mWealth(2:end,:);

end