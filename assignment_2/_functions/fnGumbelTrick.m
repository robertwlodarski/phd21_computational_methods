function LHS = fnGumbelTrick(A,B,Parameters)
    % Unpack parameters
    pZeta           = Parameters.pZeta;
    pEuler          = 0.577215;
    % Solve
    Max             = max(A/pZeta,B/pZeta);
    Inner           = exp(A / pZeta - Max) + exp(B / pZeta - Max);
    LHS             = pZeta * (pEuler + Max + log(Inner));
end