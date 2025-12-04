function [pEta]     =    fnCalibrateEta(Parameters)
    Obj             = @(x) fnGenerateMoments(setfield(Parameters, 'pEta',x)) - 0.33;
    pEta            = fsolve(Obj,1);
end