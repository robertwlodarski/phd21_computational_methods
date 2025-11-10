function [Simulations]          = fnSimulationsSolver(Parameters,Grids,SimulationsNumber)
    % 1. Unpack parameters & grids
    pT                  = Parameters.pT;
    pr                  = Parameters.pr;
    pb                  = Parameters.pb;
    pBarh               = Parameters.pBarh;
    pBeta               = Parameters.pBeta;
    pVarphi             = Parameters.pVarphi;
    vGridZ              = Grids.vGridZ;
    vGridH              = Grids.vGridH;
    vGrida              = Grids.vGrida;
    vGridAge            = Grids.vGridAge;

    % 2A. Simulate shocks drawn from U(0,1)
    rng(1997,"twister");
    mProductivityShocks = rand(pT+1,SimulationsNumber);
    mWorkingShocks      = rand(pT+1,SimulationsNumber);
    mNonWorkingShocks   = rand(pT+1,SimulationsNumber);
    mMatchingShocks     = rand(pT+1,SimulationsNumber);

    % 3. Generate matrices 
    [~,~,mW,mN]         = fnValueFunctionMatrices(Parameters,Grids);

    % 4. Prepare space
    mAssets             = zeros(size(mWorkingShocks));
    mAssetsNext         = mAssets;
    mWage               = mAssetsNext;
    mEarnings           = mWage;
    mConsumption        = mEarnings;
    mH                  = mConsumption;
    mW_opt              = mH;
    mN_opt              = mW_opt;
    mParticipation      = mN_opt;
    mMatching           = mParticipation;
    mIncome             = mMatching;
    mWorking            = mIncome;
    

    % 5. Initialise (mZ = indices; mH, mAssets = values)
    h0                  = vGridH(1);
    a0                  = 1;
    mH(1,:)             = h0;
    mAssets(1,:)        = a0;
    mZ                  = fnGenerateProductivityValues(mProductivityShocks, Grids);
    mMatching           = (mMatchingShocks < pVarphi);

    % 6. Prepare function
    F_W                 = griddedInterpolant({vGrida, vGridH, vGridZ, vGridAge, vGrida}, mW);
    F_N                 = griddedInterpolant({vGrida, vGridH, vGridZ, vGridAge, vGrida}, mN);
    fnOptions           = optimset('TolX', 1e-4, 'Display', 'off');
    
    % 7. Start simulating
    for iii = 1:1:SimulationsNumber
        for ttt = 1:1:pT+1
            % Worker states
            a                       = mAssets(ttt,iii);
            h                       = mH(ttt,iii);
            z                       = vGridZ(mZ(ttt,iii));
            w                       = fnWage(h,z,Parameters);
            ap_upper_w              = (1+pr) * a + w;
            ap_upper_n              = (1+pr) * a + pb;
            % Solve (working)
            obj_W                   = @(x) -F_W(a,h,z,vGridAge(ttt),x);
            [ap_w, W]               = fminbnd(obj_W, min(vGrida), ap_upper_w, fnOptions);
            W                       = -W;
            mW_opt(ttt,iii)         = W;
            % Solve (non-working)
            obj_N                   = @(x) -F_N(a,h,z,vGridAge(ttt),x);
            [ap_n, N]               = fminbnd(obj_N, min(vGrida), ap_upper_n, fnOptions);
            N                       = -N;
            mN_opt(ttt,iii)         = N;
            % Working probabilities 
            [dv,ds]                 = fnGumbelTrickProbabilities(W,N,Parameters);
            CondLabSupp_w           = (mWorkingShocks(ttt,iii)      <= dv);
            CondLabSupp_n           = (mNonWorkingShocks(ttt,iii)   <= ds);
            % Updating decisions and variables going into the next period
            if (ttt == 1 || mWorking(ttt-1,iii)==0) 
                % Update labour status
                mParticipation(ttt,iii)     = CondLabSupp_n;
                mWorking(ttt,iii)           = CondLabSupp_n * mMatching(ttt,iii);
                % Update asset decision
                mAssetsNext(ttt,iii)        = ap_w * mWorking(ttt,iii) + (1-mWorking(ttt,iii)) * ap_n;
                % Other things worth saving
                mWage(ttt,iii)              = w * mWorking(ttt,iii);
                mEarnings(ttt,iii)          = w * mWorking(ttt,iii) + (1-mWorking(ttt,iii)) * pb;
                mIncome(ttt,iii)            = mEarnings(ttt,iii) + (1+pr) * a;
                mConsumption(ttt,iii)       = mIncome(ttt,iii) - mAssetsNext(ttt,iii);
                % Update next skill level
                mH(ttt,iii)                 = min(h + mWorking(ttt,iii), pBarh);
            elseif mWorking(ttt-1,iii)==1 %(Working households)
                % Update labour status
                mParticipation(ttt,iii)     = CondLabSupp_w;
                mWorking(ttt,iii)           = 1;
                % Update asset decision
                mAssetsNext(ttt,iii)        = ap_w * mWorking(ttt,iii) + (1-mWorking(ttt,iii)) * ap_n;
                % Other things worth saving
                mWage(ttt,iii)              = w * mWorking(ttt,iii);
                mEarnings(ttt,iii)          = w * mWorking(ttt,iii) + (1-mWorking(ttt,iii)) * pb;
                mIncome(ttt,iii)            = mEarnings(ttt,iii) + (1+pr) * a;
                mConsumption(ttt,iii)       = mIncome(ttt,iii) - mAssetsNext(ttt,iii);
                % Update next skill level
                mH(ttt,iii)                 = min(h + mWorking(ttt,iii), pBarh);
            end 
            % Save update assets between periods
            if ttt <= pT
                mAssets(ttt+1,iii)            = mAssetsNext(ttt,iii);
            elseif ttt == pT+1
                mAssetsNext(ttt,iii)        = pBeta / (1 + pBeta) * mIncome(ttt,iii);
            end 
        end
    end

    % 8. Save the results
    Simulations.mAssets         = mAssets;
    Simulations.mAssetsNext     = mAssetsNext;
    Simulations.mWage           = mWage;
    Simulations.mEarnings       = mEarnings;
    Simulations.mIncome         = mIncome;
    Simulations.mParticipation  = mParticipation;
    Simulations.mConsumption    = mConsumption;
    Simulations.mH              = mH;
    Simulations.mWorking        = mWorking;
end