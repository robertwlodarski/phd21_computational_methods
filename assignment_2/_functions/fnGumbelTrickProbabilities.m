function [dv,ds] = fnGumbelTrickProbabilities(W,N,Parameters)
    % Unpack parameters
    pZeta           = Parameters.pZeta;
    pPhi            = Parameters.pPhi;
    pVarphi         = Parameters.pVarphi;
    % Handbaek's trick: X = max(W,N) (of them all)
    X = max(W,N) / pZeta;
    % Probability of working given working
    A               = W / pZeta;
    B               = N / pZeta;
    dv              = exp(A - X - log(exp(A-X) + exp(B-X)));

    % Probability of working given not working
    C               = (pVarphi * W + (1 - pVarphi) * N - pPhi);
    ds              = exp(C - X - log(exp(C-X) + exp(B-X)));
end