% Plotting key model statistics 

%% 1. Average values
Plotting.vAverageEmployment         = mean(Simulations.mParticipation,2);
Plotting.vAverageWealth             = mean(Simulations.mAssets,2);
Plotting.vAverageConsumption        = mean(Simulations.mConsumption,2);

%% 2. Standard plots 

figure(1);
plot(Grids.vGridAge(1:end-1),Plotting.vAverageEmployment);
grid on;
xlabel('Age');
ylabel('Proportion');