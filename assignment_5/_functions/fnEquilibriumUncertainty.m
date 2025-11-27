function [Results]          = fnEquilibriumUncertainty(Parameters, Grids,ResultsSS)

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
vGridA              = Grids.vGridA;
mTransitionA        = Grids.mTransitionA;

% Set-up elements
rng(1997);
pT                  = 20000;
pBurnIn             = 500;
pRequiredTime       = pT + pBurnIn;
iError2             = 10;
pErrorTol           = 1e-6;
vFuture             = [(2:pRequiredTime)'; pRequiredTime];
pStepSize           = 0.8;
iIterationNum       = 1;

% Set up productivity shocks
vMarginalTransition = (mTransitionA')^100 * ones(size(mTransitionA,1),1) / sum(ones(size(mTransitionA,1),1));
vCMT                = cumsum(vMarginalTransition);
vAShock             = rand(pRequiredTime,1);
vAIndex             = zeros(pRequiredTime,1);
for ttt = 1:1:pRequiredTime
    vAIndex(ttt)    = sum(vAShock(ttt) >= vCMT)+1;
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
vW                  = ones(pRequiredTime,1)*ResultsSS.W;
vN                  = ones(pRequiredTime,1)*ResultsSS.N;
vK                  = max(ones(pRequiredTime,1)*ResultsSS.K + normrnd(0,1e-8,pRequiredTime,1),1e-8);
vC                  = ones(pRequiredTime,1)*ResultsSS.C;

% Start the outer loop
tic;
while iError2 > pErrorTol

    % Update the Kp vector & "sanitise"
    vK                  = max(real(vK),1e-5);
    vN                  = max(real(vN),1e-5);
    vC                  = max(real(vC),1e-5);
    vKp                 = [vK(2:end);vK(1)];
    
    %% 3. Backward solution
    % [No time loop; vectorisation used]
    iEV1                = 0.0;

    % Run the process for every A (vectorised for all Ks)
    for iApIndex        = 1:1:size(vGridA,1) 

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
        iKLow                                   = sum(repmat(iCandidate',length(vKp),1)< vKp,2);
        iKLow(iKLow < 1)                        = 1;
        iKLow(iKLow >= length(iIndex))          = length(iIndex)-1;
        iKHigh                                  = iKLow + 1;
        iWeightLow                              = (iCandidate(iKHigh) - vKp) ./ (iCandidate(iKHigh) - iCandidate(iKLow));
        iWeightLow(iWeightLow < 0)              = 0;
        iWeightLow(iWeightLow > 1)              = 1;

        % D. Extrapolation for wages, N, and adjustment cost derivative
        iWLow                                   = (1 - pAlpha) * Ap * vK(iKLow).^(pAlpha) .* vN(iKLow).^(-pAlpha);
        iWHigh                                  = (1 - pAlpha) * Ap * vK(iKHigh).^(pAlpha) .* vN(iKHigh).^(-pAlpha);
        iNLow                                   = (vC(iKLow).^(-pSigma) .* iWLow / pEta).^(pChi);
        iNHigh                                  = (vC(iKHigh).^(-pSigma) .* iWHigh / pEta).^(pChi);
        iRLow                                   = pAlpha * Ap * vK(iKLow).^(pAlpha-1).*iNLow.^(1-pAlpha)-pDelta;
        iRHigh                                  = pAlpha * Ap * vK(iKHigh).^(pAlpha-1).*iNHigh.^(1-pAlpha)-pDelta;
        iPsi1Low                                = -pMu * (vKp(iKLow)-vK(iKLow)) ./ (vK(iKLow)) .* vKp(iKLow)./vK(iKLow)+ pMu / 2 * ((vKp(iKLow) - vK(iKLow))./vK(iKLow)).^2;
        iPsi1High                               = -pMu * (vKp(iKHigh)-vK(iKHigh)) ./ (vK(iKHigh)) .* vKp(iKHigh)./vK(iKHigh) + pMu / 2 * ((vKp(iKHigh) - vK(iKHigh))./vK(iKHigh)).^2;

        % E. Expected derivative of value
        iV1Low                                  = vC(iKLow).^(-pSigma) .* (1 + iRLow - iPsi1Low);
        iV1High                                 = vC(iKHigh).^(-pSigma) .* (1 + iRHigh - iPsi1High);
        iEV1                                    = iEV1 + pBeta * (vAp ~= Ap) .* mTransitionA(vAIndex,iApIndex) .* (iWeightLow .* iV1Low + (1-iWeightLow) .* iV1High);
    end

    % F. Add the realised path
    iWFuture            = (1 - pAlpha) * vAp .* vKp.^(pAlpha) .* vN(vFuture).^(-pAlpha);
    iNFuture            = (vC(vFuture).^(-pSigma) .* iWFuture / pEta).^(pChi);
    iRFuture            = pAlpha * vAp .* vKp.^(pAlpha-1).*iNFuture.^(1-pAlpha)-pDelta;
    iPsiFuture          = -pMu * (vKp(vFuture)-vKp) ./ (vKp).* vKp(vFuture)./vKp + pMu / 2 * ((vKp(vFuture) - vKp)./vKp).^2;
    iVFuture            = vC(vFuture).^(-pSigma) .* (1 + iRFuture - iPsiFuture);
    iEV1                = iEV1 + pBeta * vRealisedP .* iVFuture;

    % G. Update the allocations
    iDenominator        = 1 + pMu * (vKp - vK) ./ vK;
    vCTemp              = (iEV1 ./ iDenominator).^(-1 / pSigma);
    vN                  = (vCTemp.^(-pSigma) .* vW / pEta).^(pChi);
    vI                  =  vA .* vK .^(pAlpha) .* vN .^(1 - pAlpha) -  pMu / 2 * ((vKp - vK) ./ vK).^2 -vCTemp;

    %% 4. Iterate forward
    vKPast              = [vK(end);vK(1:end-1)];
    vIPast              = [vI(end);vI(1:end-1)];
    vKNew               = vKPast * (1 - pDelta) + vIPast;
    vKpNew              = [vKNew(2:end);vKNew(1)];
    vINew               = vKpNew - (1 - pDelta) * vKNew;
    vCNew               = vA .* vKNew .^(pAlpha) .* vN .^(1 - pAlpha) -  pMu / 2 * ((vKpNew - vKNew) ./ vKNew).^2 - vINew;

    % Compute errror
    iError2             = mean([vC-vCNew;vK-vKNew].^2,"all");
    iIterationNum       = iIterationNum + 1;

    % Update guesses
    vC                  = pStepSize * vC + (1 - pStepSize) * vCNew;
    vK                  = pStepSize * vK + (1 - pStepSize) * vKNew;

    % Update other used values
    vW                  = ((1 - pAlpha) * vA .* vK.^(pAlpha) .* (vC .^(-pSigma) / pEta).^(- pAlpha * pChi)).^(1 / (1 + pAlpha * pChi));
    vN                  = (vC.^(-pSigma) .* vW / pEta).^(pChi);

    %% 5. Optional reports
    if (floor((iIterationNum-1)/25) == (iIterationNum-1)/25)
        fprintf('Convergence: \n');
        fprintf('Iteration: %.0f \n', iIterationNum);
        fprintf('Error:     %.8f \n', iError2);
    else 
    end

% End the outer loops
end 
toc;

%% Save the results
Results.vK              = vK;

end 