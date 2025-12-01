function [Hours]  = fnGenerateMoments(Parameters)
% Run the model
Results            = fnSolveSSEquilibrium(Parameters);
Hours              = Results.N;
end