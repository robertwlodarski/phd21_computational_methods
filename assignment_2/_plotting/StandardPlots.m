% Plotting key model statistics 

%% 1. Average values
Plotting.vAverageEmployment         = mean(Simulations.mParticipation,2);
Plotting.vAverageWealth             = mean(Simulations.mAssets,2);
Plotting.vAverageConsumption        = mean(Simulations.mConsumption,2);

%% 2. Standard plots 

% Employment
figure(1);
plot(Grids.vGridAge(1:end-1),Plotting.vAverageEmployment,'LineWidth',2.5);
grid on;
xlabel('Age');
ylabel('Proportion');
saveas(gcf,'_figures/AverageEmployment.png');

% Wealth
figure(2);
plot(Grids.vGridAge,Plotting.vAverageWealth,'LineWidth',2.5);
grid on;
xlabel('Age');
ylabel('Proportion');
saveas(gcf,'_figures/AverageWealth.png');

% Consumption
figure(3);
plot(Grids.vGridAge(1:end-1),Plotting.vAverageConsumption,'LineWidth',2.5);
grid on;
xlabel('Age');
ylabel('Proportion');
saveas(gcf,'_figures/AverageConsumption.png');