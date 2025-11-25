function [Results]          = fnEquilibriumUncertainty(Parameters, Grids,ResultsSS)

%% 1. Prepare key elements

% From parameters
pSigma              = Parameters.pSigma;
pEta                = Parameters.pEta;
pChi                = Parameters.pChi;
pBeta               = Parameters.pBeta;

% From grids
vGridA              = Grids.vGridA;
mTransitionA        = Grids.mTransitionA;

% Set-up
rng(1997);
pT                  = 5000;
pBurnIn             = 500;
pRequiredTime       = pT + pBurnIn;

%% 2. Initial guesses

% Run the steady state version
vMarginalTransition = (mTransitionA')^100 * ones(size(mTransitionA,1),1) / sum(ones(size(mTransitionA,1),1));
pEA                 = vMarginalTransition' * vGridA;
pEV                 = (log(ResultsSS.C) - pEta * (ResultsSS.N)^(1 + 1/pChi)/(1 + 1/pChi)) / (1 - pBeta);

% Initial guesses
vR                  = ones(pRequiredTime,1)*ResultsSS.R;
vN                  = ones(pRequiredTime,1)*ResultsSS.N;
vK                  = ones(pRequiredTime,1)*ResultsSS.K + normrnd(0,1e-5,pRequiredTime,1);
vV                  = ones(pRequiredTime,1)*pEV;

end 