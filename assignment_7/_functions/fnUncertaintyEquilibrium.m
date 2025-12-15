function [Results]          = fnUncertaintyEquilibrium(Parameters, Grids,Results0)

%% 1. Prepare key elements

% From parameters
pSigma              = Parameters.pSigma;
pEta                = Parameters.pEta;
pChi                = Parameters.pChi;
pBeta               = Parameters.pBeta;
pAlpha              = Parameters.pAlpha;
pDelta              = Parameters.pDelta;
pMu                 = Parameters.pMu;

% From grids
vGridA              = Grids.vGridTFP;
mTransitionA        = Grids.mTransitionTFP;
vGridZ              = Grids.vGridZ;
vGridA2             = Grids.vGridA2;

% Set-up elements
rng(1997);
pT                  = 10000;
pBurnIn             = 500;
pRequiredTime       = pT + pBurnIn;
iError2             = 10;
pErrorTol           = 1e-8;
vFuture             = [(2:pRequiredTime)'; pRequiredTime];
pStepSize           = 0.95;
iIterationNum       = 1;

% Set up productivity shocks
vAShock             = rand(pRequiredTime,1);
vAIndex             = zeros(pRequiredTime,1);
vAIndex(1)          = 4;
for ttt = 2:1:pRequiredTime
    vAIndex(ttt)    = sum(vAShock(ttt) >= cumsum(squeeze(mTransitionA(vAIndex(ttt-1),:))))+1;
end 
vApIndex            = [vAIndex(2:end);vAIndex(1)];
vA                  = vGridA(vAIndex);
vAp                 = vGridA(vApIndex);
vRealisedP          = zeros(pRequiredTime,1);
for ttt = 1:1:pRequiredTime
    vRealisedP(ttt) = mTransitionA(vAIndex(ttt),vApIndex(ttt));
end

%% 2. Initial guesses

% Initial guesses 
% Aggregate
vL                  = ones(pRequiredTime,1)*Results0.vLabourSupply;
vK                  = max(ones(pRequiredTime,1)*Results0.vCapitalOpt + normrnd(0,1e-8,pRequiredTime,1),1e-8);

% Policy function (C, Ap, and N)
mPolicyC                    = repmat(Results0.mC,1,1,pRequiredTime);
mPolicyN                    = repmat(Results0.mN,1,1,pRequiredTime);
mPolicyA2                   = repmat(Results0.mPolicyWealthNext,1,1,pRequiredTime);

%% 3. Start running the loop

% Start the outer loop
tic;
while iError2 > pErrorTol

    %% 4. Compute the price guesses
    vKtoL               = vK ./ vL;
    vInterest           = - pDelta + pAlpha * vA .* vKtoL.^(pAlpha-1);
    vWage               = (1-pAlpha) * vA .* vKtoL.^(pAlpha);
    vKp                 = [vK(2:end);vK(1)];
    vLp                 = [vL(2:end);vL(1)];
    vKptoLp             = vKp./vLp;
    vInterestFuture     = - pDelta + pAlpha * vAp .* vKptoLp.^(pAlpha-1);

    %% 5. Backward simulation
    % Set the space
    mEV1                = zeros(size(vGridA2,1),size(vGridZ,1),pRequiredTime);

    % Open the TFP loop
    for iApIndex                = 1:1:size(vGridA,1)

        %% 5.1. Find the "neighbours"
        % A. Look for candidates
        Ap                                      = vGridA(iApIndex);
        iCandidateLoc                           = find(vA == Ap);

        % B. Clean "forbidden" candidates & sort
        iCandidateLoc(iCandidateLoc<pBurnIn)    = [];
        iCandidateLoc(iCandidateLoc>pT)         = [];
        iCandidate                              = vK(iCandidateLoc);
        [iCandidate, iIndex]                    = sort(iCandidate);
        iCandidateLoc                           = iCandidateLoc(iIndex);

        % C. Extrapolation setting
        iEdges                                  = [-Inf; iCandidate; Inf];
        iBinIndices                             = discretize(vKp, iEdges);
        iKLow                                   = iBinIndices -1;
        iKLow(iKLow < 1)                        = 1;
        iKLow(iKLow >= length(iIndex))          = length(iIndex)-1;
        iKHigh                                  = iKLow + 1;
        iIndLow                                 = iCandidateLoc(iKLow);
        iIndHigh                                = iCandidateLoc(iKHigh);
        iWeightLow                              = (vK(iIndLow) - vKp) ./ (vK(iIndHigh) - vK(iIndLow));
        iWeightLow(iWeightLow < 0)              = 0;
        iWeightLow(iWeightLow > 1)              = 1;

        % D. Extrapolate interest rates
        iRLow                                   = pAlpha * Ap * vKtoL(iIndLow).^(pAlpha-1)-pDelta;
        iRHigh                                  = pAlpha * Ap * vKtoL(iIndHigh).^(pAlpha-1)-pDelta;

        for iZIndex             = 1:1:size(vGridZ,1)
            
            %% 5.2. Obtain the EV components
            % A. Policy functions
            iCLow               = squeeze(mPolicyC(:,iZIndex,iIndLow));
            iCHigh              = squeeze(mPolicyC(:,iZIndex,iIndHigh));
            iapLow              = squeeze(mPolicyA2(:,iZIndex,iIndLow));
            iapHigh             = squeeze(mPolicyA2(:,iZIndex,iIndHigh));

            % B. Generate labour responses and adjustment costs
            iPsi2Low            = pMu / 2 * (1 - (iapLow.^2./vGridA2.^2));
            iPsi2High           = pMu / 2 * (1 - (iapHigh.^2./vGridA2.^2));

            % C. Generate the lower and upper expected value function
            iEVLow              = (1 + repmat(iRLow',size(vGridA2,1),1)-iPsi2Low) ./ iCLow;
            iEVHigh             = (1 + repmat(iRHigh',size(vGridA2,1),1)-iPsi2High) ./ iCHigh;
            iEVInterpolated     = repmat(iWeightLow',size(vGridA2,1),1) .* iEVLow + (1 - repmat(iWeightLow',size(vGridA2,1),1)) .* iEVHigh;
            
            % D. Obtain the value 
            mEV1(:,iZIndex,:)   = squeeze(mEV1(:,iZIndex,:)) + pBeta * repmat(((vAp ~= Ap) .* mTransitionA(vAIndex,iApIndex))',size(vGridA2,1),1) .* iEVInterpolated;

        % Close the skill shock loop
        end
        % Close the TFP loop
    end

    % E. Add the realised path
    iCFuture                    = mPolicyC(:,:,[2:end,1]);
    iPsi2Future                 = pMu / 2 * (1 - (mPolicyA2.^2 ./ vGridA2.^2));
    iPsi2Future                 = iPsi2Future(:,:,[2:end,1]);
    iVFuture                    = (1+ repmat(permute(vInterestFuture,[2,3,1]),size(vGridA2,1),size(vGridZ,1),1) - iPsi2Future)./iCFuture;
    mEV1                        = mEV1 + pBeta *  repmat(permute(vRealisedP,[2,3,1]),size(vGridA2,1),size(vGridZ,1),1) .* iVFuture;

    %% 5.3. Update the allocations
    % A. New, temporary consumption level
    iPsi1                       = pMu * ((mPolicyA2 ./ vGridA2-1));
    mCTempUnc                   = (1 + iPsi1) ./ mEV1;
    mPolicyCUpdated             = max(1e-3,mCTempUnc);

    % B. Updated labour and asset policies 
    mPolicyNUpdated             =  (repmat(permute(vWage,[2,3,1]),size(vGridA2,1),size(vGridZ,1),1) .* repmat(vGridZ',size(vGridA2,1),1,pRequiredTime) ./ (pEta * mPolicyCUpdated));
    mPsi                        = pMu / 2 * (mPolicyA2.^2 - 2 * mPolicyA2 .* repmat(vGridA2,1,size(vGridZ,1),pRequiredTime) + repmat(vGridA2,1,size(vGridZ,1),pRequiredTime).^2) ./repmat(vGridA2,1,size(vGridZ,1),pRequiredTime);
    mPolicyA2Updated            = (1+repmat(permute(vInterest,[2,3,1]),size(vGridA2,1),size(vGridZ,1),1)).* repmat(vGridA2,1,size(vGridZ,1),pRequiredTime)...                   %Interest income
                                    +mPolicyNUpdated.*repmat(permute(vWage,[2,3,1]),size(vGridA2,1),size(vGridZ,1),1) .* repmat(vGridZ',size(vGridA2,1),1,pRequiredTime)...     %Labour income
                                    -mPolicyCUpdated-mPsi;
    

    %% 6. Non-stochastic simulation

% End the outer loops
end 
toc;

%% Save the results


end 