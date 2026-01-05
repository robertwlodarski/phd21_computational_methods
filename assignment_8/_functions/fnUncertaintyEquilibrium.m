function Results        = fnUncertaintyEquilibrium(Parameters,Grids, ResultSS)

%% 1. Prepare the setting

% Parameters
pEta                = Parameters.pEta;


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


%% 3. Start running the loop
tic;
while iError > pErrorTol
    fprintf('Iteration:             %.0f \n', iIterationNum);

    % Implied elements calculation
    vW                          = pEta * vC;
    vKp                         = [vK(2:end);vK(1)];

    %% 4. Backward simulation begins
    
    % 4.1. Initialise the expected marginal value
    mEJ1                        = zeros(size(vGridK,1),pRequiredTime);

    % 4.2. Open the TFP loop
    for iApIndex = 1:1:size(vGridA,1)
        
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

        % Close the TFP loop
    end 

    iIterationNum               = iIterationNum + 1;
    % Close the outer loop
end 
toc;