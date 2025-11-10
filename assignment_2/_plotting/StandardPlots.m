% Plotting key model statistics 

%% 1. Average values
Plotting.vAverageEmployment         = mean(Simulations.mWorking,2);
Plotting.vAverageWealth             = mean(Simulations.mAssets,2);
Plotting.vAverageConsumption        = mean(Simulations.mConsumption,2);

%% 2. Standard plots 

% Employment
figure(1);
plot(24+Grids.vGridAge,Plotting.vAverageEmployment,'LineWidth',2.5);
grid on;
xlabel('Age');
ylabel('Proportion');
saveas(gcf,'_figures/AverageEmployment.png');

% Wealth
figure(2);
plot(24+Grids.vGridAge,Plotting.vAverageWealth,'LineWidth',2.5);
grid on;
xlabel('Age');
ylabel('');
saveas(gcf,'_figures/AverageWealth.png');

% Consumption
figure(3);
plot(24+Grids.vGridAge,Plotting.vAverageConsumption,'LineWidth',2.5);
grid on;
xlabel('Age');
ylabel('');
saveas(gcf,'_figures/AverageConsumption.png');