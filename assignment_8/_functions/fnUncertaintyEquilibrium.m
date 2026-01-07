function Results        = fnUncertaintyEquilibrium(Parameters,Grids, ResultSS)

%% 1. Prepare the setting
% Parameters
pEta                = Parameters.pEta;
pAlpha              = Parameters.pAlpha;
pDelta              = Parameters.pDelta;
pGamma              = Parameters.pGamma;
pBeta               = Parameters.pBeta;
pPhi                = Parameters.pPhi;

% Grids
vGridA              = Grids.vGridTFP;
vGridK              = Grids.vGridK;
mTransitionA        = Grids.mTransitionTFP;

% Set-up elements
rng(1997);
pT                  = 10000;
pBurnIn             = 500;
pRequiredTime       = pT + pBurnIn;
iError              = 10;
pErrorTol           = 1e-4;
pStepSize           = 0.99;
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

% Guess capital (with minor pertubation) and consumption
vK                  = max(ResultSS.K * ones(pRequiredTime,1) + normrnd(0,1e-8,pRequiredTime,1),1e-8);
vC                  = ResultSS.C * ones(pRequiredTime,1);
vIMin               = pPhi * pDelta * ResultSS.K;

% Steady state value function
vInnerTerm          = pBeta / ResultSS.C * repmat(pAlpha * ResultSS.Y / ResultSS.K + 1 - pDelta,pRequiredTime,1);

%% 3. Start running the loop
tic;
while iError > pErrorTol
    fprintf('Iteration:             %.0f \n', iIterationNum);

    % Implied elements calculation
    vW                          = pEta * vC;
    vKp                         = [vK(2:end);vK(1)];
    vN                          = ((pGamma * vA .* vK.^(pAlpha))./vW).^(1 / (1-pGamma));

    %% 4. Backward simulation begins
    
    % 4.1. Initialise the expected marginal value
    vEJ1                        = 0.0;

    % 4.2. Open the TFP loop
    for iApIndex = 1:1:size(vGridA,1)
        
        %% 5. Find the "neighbours"
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

        % D. Get MPK values for unconstrained solution
        iInnerTermLow                           = vInnerTerm(iIndLow); 
        iInnerTermHigh                          = vInnerTerm(iIndHigh); 
        iInnerInterpolated                      = iWeightLow .* iInnerTermLow + (1 - iWeightLow) .* iInnerTermHigh;
        vEJ1                                    = vEJ1 +  mTransitionA(vAIndex,iApIndex) .* iInnerInterpolated;

        % Close the TFP loop
    end 
    
    %% 6. Update the allocations

    % 6.1. Unconstrained allocation
    vY                          = vA .* vK.^(pAlpha) .* vN.^(pGamma);
    vCUnc                       = 1 ./ vEJ1;
    vKpUnc                      = vY + (1 - pDelta) .* vK - vCUnc;

    % 6.2. Add constraints and obtain implied K and C
    vKpMin                      = vIMin + (1-pDelta)*vK;
    vKpImplied                  = max(vKpMin,vKpUnc);
    vCImplied                   =  vY + (1 - pDelta) .* vK - vKpImplied;
    vWImplied                   = pEta * vCImplied;
    vNImplied                   = (pGamma * vA .* vK.^(pAlpha) ./ vWImplied).^(1 / (1 - pGamma));
    vYImplied                   = vA .* vK.^(pAlpha) .* vNImplied.^(pGamma);

    % 6.3. Compute the inner term for the next iteration
    vCpImplied                  = [vCImplied(2:end);vCImplied(1)];
    vShadowPrice                = min(vCImplied .* vEJ1,1.0);
    vInnerTerm                  = pBeta ./ vCImplied .* (pAlpha * vYImplied ./ vK + (1-pDelta) .* vShadowPrice);

    %% 7. Finalise the loop

    % Adjust vK
    vKImplied                   = [vKpImplied(end);vKpImplied(1:end-1)];

    % Compute errror
    iError                      = mean([vC-vCImplied;vK-vKImplied].^2,"all");
    iErrorK                     = mean((vK-vKImplied).^2,"all");
    iErrorC                     = mean((vC-vCImplied).^2,"all");

    % Update the logic
    vK                          = pStepSize * vK + (1 - pStepSize)*vKImplied;
    vC                          = pStepSize * vC + (1 - pStepSize)*vCImplied;
    
     %% 8. Optional reports
     % Print out
     if (floor((iIterationNum-1)/20) == (iIterationNum-1)/20)
         fprintf('====================================== \n');
         fprintf('====================================== \n');
         fprintf('Convergence: \n');
         fprintf('Iteration:         %.0f \n', iIterationNum);
         fprintf('Error:             %.8f \n', iError);
         fprintf('Capital error:     %.8f \n', iErrorK);
         fprintf('Consumption error: %.8f \n', iErrorC);
     else
     end
    iIterationNum               = iIterationNum + 1;
    % Close the outer loop
end 
toc;

%% 9. Save the results
Results.vK                      = vK;
Results.vC                      = vC;
Results.vA                      = vA;
Results.vY                      = vYImplied;
Results.vShadowPrice            = vShadowPrice;