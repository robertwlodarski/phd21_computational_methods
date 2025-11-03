function [Simulations]          = fnSimulationsSolver(Parameters,Grids,SimulationsNumber)
    % 1. Unpack parameters of interest
    pT                          = Parameters.pT;

    % 2. Simulate shocks drawn from U(0,1)
    rng(1997,"twister");
    vProductivityShocks         = rand(pT,SimulationsNumber);
    vWorkingShocks              = rand(pT,SimulationsNumber);
    vNonWorkingShocks           = rand(pT,SimulationsNumber);

    % 3. Generate matrices 
    [mV,mS,mW,mN]      = fnValueFunctionMatrices(Parameters,Grids);

    % 4. Prepare space
    vH                  = zeros(pT, SimulationsNumber);
    vAssets             = vH;
    vWorkingInd         = vAssets;
    vWage               = vWorkingInd;
    vConsumption        = vWage;

    % 5. Initialise values 
    h0                  = 0;
    a0                  = 1;
end