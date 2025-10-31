function LHS = fnGumbelTrick(A,B,Parameters)
    % Unpack parameters
    pZeta           = Parameters.pZeta;
    pEuler          = 0.577215;
    % Solve
    Inner           = exp(A / pZeta) + exp(B / pZeta);
    LHS             = pZeta * (pEuler + log(Inner));
end