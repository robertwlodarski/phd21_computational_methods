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
%xline(50,'--','LineWidth',1.5,'red');
xlabel('$t$','interpreter','latex','FontSize',14);
ylabel('$Ea$','interpreter','latex','FontSize',14);
legend('Simulations','Histogram','Location','best');
saveas(gcf,'_figures/SimulationsWealth.png');