%% 1. Average wealth plot

% Data
Plots.SimAvgWealth      = mean(ResultsSimulation.mAssets,2);
Plots.HisAvgWalth       = Grids.vGridA2' * sum(Results.mDistribution,2);
Plots.HisAvgWalth       = repmat(Plots.HisAvgWalth,Parameters.pT,1);

% Figure
figure(1);
plot(1:1:Parameters.pT, Plots.SimAvgWealth, 'LineWidth', 2.5);
hold on;
grid on;
plot(1:1:Parameters.pT, Plots.HisAvgWalth, 'LineWidth', 2.5);
xline(50,'--','Color','r','LineWidth',1.5);
xlabel('$t$','interpreter','latex','FontSize',14);
ylabel('$Ea$','interpreter','latex','FontSize',14);
legend('Simulations','Histogram','Location','best');
saveas(gcf,'_figures/SimulationsWealth.png');

%% 2. 95th perctentile plot

Plots.Sim95Wealth       = prctile(ResultsSimulation.mAssets,95,2);

% Figure
figure(2);
plot(1:1:Parameters.pT, Plots.Sim95Wealth, 'LineWidth', 2.5);

