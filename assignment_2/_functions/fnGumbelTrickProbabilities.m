function [dv,ds] = fnGumbelTrickProbabilities(W,N,Parameters)
    % Unpack parameters
    pZeta           = Parameters.pZeta;
    pPhi            = Parameters.pPhi;
    pVarphi         = Parameters.pVarphi;
    % Tricking Matlab to behave well
    X1              = min(W,N);
    X2              = min((pVarphi * W + (1 - pVarphi) * N - pPhi),N);
    % Probability of working given working
    Den1            = exp( (W-X1) / pZeta ) + exp( (N-X1) / pZeta);
    dv              = exp( (W-X1) / pZeta) / Den1;
    % Probability of working given not working
    Inner           = (pVarphi * W + (1 - pVarphi) * N - pPhi);
    Den2            = exp( (Inner-X2) / pZeta) + exp( (N-X2) / pZeta);
    ds              = exp( (Inner-X2) / pZeta ) / Den2;
end