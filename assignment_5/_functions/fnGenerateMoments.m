function [Hours]  = fnGenerateMoments(Parameters, Grids)
% Run the model
Results            = fnSteadyStateSolverPFI(1.0,Parameters,Grids);
Hours              = Results.N;
end