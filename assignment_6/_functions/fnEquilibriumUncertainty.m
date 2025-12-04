function [Results]          = fnEquilibriumUncertainty(Parameters, Grids,ResultsSS)

%% 1. Prepare key elements

% From parameters
pSigma              = Parameters.pSigma;
pEta                = Parameters.pEta;
pChi                = Parameters.pChi;
pBeta               = Parameters.pBeta;
pAlpha              = Parameters.pAlpha;
pDelta              = Parameters.pDelta;
pPsi                = Parameters.pPsi;
pKss                = ResultsSS.K;

% From grids
vGridA              = Grids.vGridA;
mTransitionA        = Grids.mTransitionA;

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
vC                  = ones(pRequiredTime,1)*ResultsSS.C;
vK                  = max(ones(pRequiredTime,1)*ResultsSS.K + normrnd(0,1e-8,pRequiredTime,1),1e-8);

% Start the outer loop
tic;
while iError2 > pErrorTol

    % Update the Kp vector
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

        % D. Extrapolation for wages, N, and adjustment cost derivative
        iNLow                                   = (((1-pAlpha) * Ap .* vK(iIndLow).^(pAlpha) .* vC(iIndLow).^(-pSigma)) / pEta).^(pChi / (1+ pAlpha* pChi));
        iNHigh                                  = (((1-pAlpha) * Ap .* vK(iIndHigh).^(pAlpha) .* vC(iIndHigh).^(-pSigma)) / pEta).^(pChi / (1+ pAlpha* pChi));
        iRLow                                   = pAlpha * Ap * vK(iIndLow).^(pAlpha-1).*iNLow.^(1-pAlpha)-pDelta;
        iRHigh                                  = pAlpha * Ap * vK(iIndHigh).^(pAlpha-1).*iNHigh.^(1-pAlpha)-pDelta;

        % E. Expected derivative of value
        iV1Low                                  = vC(iIndLow).^(-pSigma) .* (1 + iRLow);
        iV1High                                 = vC(iIndHigh).^(-pSigma) .* (1 + iRHigh);
        iEV1                                    = iEV1 + pBeta * (vAp ~= Ap) .* mTransitionA(vAIndex,iApIndex) .* (iWeightLow .* iV1Low + (1-iWeightLow) .* iV1High);
    end

    % F. Add the realised path
    iNFuture            = (((1-pAlpha) * vAp .* vKp.^(pAlpha) .* vC(vFuture).^(-pSigma)) / pEta).^(pChi / (1+ pAlpha* pChi));
    iRFuture            = pAlpha * vAp .* vKp.^(pAlpha-1).*iNFuture.^(1-pAlpha)-pDelta;
    iVFuture            = vC(vFuture).^(-pSigma) .* (1 + iRFuture);
    iEV1                = iEV1 + pBeta * vRealisedP .* iVFuture;

    % G1. Update the allocations [unconstrained]
    vCTempUnc           = (iEV1).^(-1 / pSigma);
    vNUnc               = (((1-pAlpha) * vA .* vK.^(pAlpha) .* vCTempUnc.^(-pSigma)) / pEta).^(pChi / (1+ pAlpha* pChi));
    vYUnc               = vA .* vK.^(pAlpha) .* vNUnc.^(1 - pAlpha);
    vIUnc               = vYUnc - vCTempUnc;

    % G2. Update the allocations [constrained]
    vConsIndicator      = (vIUnc < pPsi * pDelta * pKss);
    vI                  = vIUnc .* (~ vConsIndicator) + pPsi * pDelta * pKss * (vConsIndicator);
    vCTemp              = vCTempUnc .* (~ vConsIndicator) + (vYUnc -  pPsi * pDelta * pKss) .*(vConsIndicator);
    vN                  = (((1-pAlpha) * vA .* vK.^(pAlpha) .* vCTemp.^(-pSigma)) / pEta).^(pChi / (1+ pAlpha* pChi));

    %% 4. Iterate forward
    vKPast              = [vK(end);vK(1:end-1)];
    vIPast              = [vI(end);vI(1:end-1)];
    vKNew               = vKPast * (1 - pDelta) + vIPast;
    vKpNew              = [vKNew(2:end);vKNew(1)];
    vINew               = vKpNew - (1 - pDelta) * vKNew;
    vCNew               = vA .* vKNew .^(pAlpha) .* vN .^(1 - pAlpha) - vINew;

    % Compute errror
    iError2             = mean([vC-vCNew;vK-vKNew].^2,"all");
    iErrorK             = mean((vK-vKNew).^2,"all");
    iErrorC             = mean((vC-vCNew).^2,"all");
    iIterationNum       = iIterationNum + 1;

    % Update guesses
    vC                  = pStepSize * vC + (1 - pStepSize) * vCNew;
    vK                  = pStepSize * vK + (1 - pStepSize) * vKNew;

    %% 5. Optional reports
    % Print out
    if (floor((iIterationNum-1)/50) == (iIterationNum-1)/50)
        fprintf('====================================== \n');
        fprintf('====================================== \n');
        fprintf('Convergence: \n');
        fprintf('Iteration:         %.0f \n', iIterationNum);
        fprintf('Error:             %.8f \n', iError2);
        fprintf('Capital error:     %.8f \n', iErrorK);
        fprintf('Consumption error: %.8f \n', iErrorC);
        fprintf('------------------------------------- \n');
        fprintf('Mean consumption:  %.2f \n', mean(vC));
        fprintf('Mean capital:      %.2f \n', mean(vK));
        fprintf('Mean employment:   %.2f \n', mean(vN));
    else
    end

    % Print out plots (more rarely than the report)
    if (floor((iIterationNum-1)/200) == (iIterationNum-1)/200)
        subplot(1,2,1);
        plot(pBurnIn:pRequiredTime-pBurnIn,vK(pBurnIn:pRequiredTime-pBurnIn)-vKNew(pBurnIn:pRequiredTime-pBurnIn),'LineWidth',1,'Color','r');
        grid on;
        xlim([pBurnIn,pRequiredTime-pBurnIn]);
        yline(0,'LineStyle','--','Color','black','LineWidth',1.5);
        xlabel('Time')
        ylabel('Prediced - realised K')

        subplot(1,2,2);
        plot(pBurnIn:pRequiredTime-pBurnIn,vC(pBurnIn:pRequiredTime-pBurnIn)-vCNew(pBurnIn:pRequiredTime-pBurnIn),'LineWidth',1,'Color','r');
        grid on;
        xlim([pBurnIn,pRequiredTime-pBurnIn]);
        yline(0,'LineStyle','--','Color','black','LineWidth',1.5);
        xlabel('Time')
        ylabel('Prediced - realised C')
        drawnow;
    else
    end
% End the outer loops
end 
toc;

%% Save the results
Results.vK              = vK(pBurnIn:pRequiredTime-pBurnIn);
Results.vA              = vA(pBurnIn:pRequiredTime-pBurnIn);
Results.vC              = vC(pBurnIn:pRequiredTime-pBurnIn);
Results.vI              = vI(pBurnIn:pRequiredTime-pBurnIn);

end 
