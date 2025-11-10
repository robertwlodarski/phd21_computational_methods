% Value functions decompositions

%% 1. Average values
% Data
Decompositions.mWorkingW                            = Simulations.mWopt;
Decompositions.mWorkingW(Simulations.mWorking~=1)   = NaN;

% Decompositions
Decompositions.vAverageW                  = mean(Simulations.mWopt,2);
Decompositions.vAverageN                  = mean(Simulations.mNopt,2);
Decompositions.mGumbelProbs               = zeros(Parameters.pT,2);
for iii = size(Decompositions.mGumbelProbs,1)
   Decompositions.mGumbelProbs(iii,:) = fnGumbelTrickProbabilities(Decompositions.vAverageW(iii),Decompositions.vAverageN(iii),Parameters); 
end




Decompositions.vWorkingConsumption        = log(mean(Simulations.mConsumption(Simulations.mWorking==1),2));
Decompositions.vNonWorkingConsumption     = log(mean(Simulations.mConsumption(Simulations.mWorking~=1),2));