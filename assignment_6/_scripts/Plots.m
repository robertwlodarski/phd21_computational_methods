%% 1. Comparing version with and without the ocassionally binding constraints
Parameters.pPsi         = 0.975;

% Investment
figure(1)
plot(1:1:500,Results.vI(1:1:500),'Color','blue','LineWidth',2.5);
hold on;
grid on;
plot(1:1:500,ResultsNC.vI(1:1:500),'Color','red','LineWidth',2.5);
yline(ResultsSS.K * Parameters.pPsi * Parameters.pDelta,'LineStyle','--','Color','black');
legend({'Constraint binding','No constraint'},'Location','best');
hold off;
saveas(gcf,'_figures/InvestmentComparison.png');

% Consumption
figure(2)
plot(1:1:500,Results.vC(1:1:500),'Color','blue','LineWidth',2.5);
hold on;
grid on;
plot(1:1:500,ResultsNC.vC(1:1:500),'Color','red','LineWidth',2.5);
yline(ResultsSS.C,'LineStyle','--','Color','black');
legend({'Constraint binding','No constraint'},'Location','best');
saveas(gcf,'_figures/ConsumptionComparison.png');
