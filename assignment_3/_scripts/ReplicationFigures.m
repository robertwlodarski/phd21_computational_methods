% Replication of figures 2a and 2b

%% 2a. Interest rate vs. per capita assets

% Bounds for interest
pBeta               = Parameters.pBeta;
pUpperInterst       = 1/ pBeta - 1;
vInterests          = linspace(-0.01,0.035,50);

% Linearised function
InterestAssets      = @(x) arrayfun(@(x_scalar) fnInterestAssets(x_scalar,Grids,Parameters), x);
vAssets             = InterestAssets(vInterests);

% Plot 
figure(1);
plot(vAssets,vInterests,'LineWidth',2.5);
yline(pUpperInterst, '--', '\frac{1-\beta}{\beta}','interpreter','latex','LineWidth',1.5,'FontSize',14,'LabelVerticalAlignment','bottom');
yline(0, 'r', 'LineWidth',1.5);
grid on;
xlabel('$Ea$','interpreter','latex','FontSize',14);
ylabel('$r$','interpreter','latex','FontSize',14);
saveas(gcf,'_figures/Figure2aReplication.png');

%% 2b. Equilibrium

% Parameters
pDelta              = Parameters.pDelta;
pAlpha              = Parameters.pAlpha;
pA                  = Parameters.pA;
vL                  = Results.vLabourSupply;

% Get capital
Capital             = @(x) arrayfun(@(r) ((r + pDelta) / (pAlpha*pA))^(1/(pAlpha-1))*vL ,x);
vCapital            = Capital(vInterests);

% Plot
figure(2);
plot(vInterests, vCapital, 'LineWidth', 2.5);
hold on;
plot(vInterests, vAssets, 'LineWidth', 2.5)
grid on;
xlabel('$r$','interpreter','latex','FontSize',14);
ylabel('$K$','interpreter','latex','FontSize',14);
legend('Capital demand','Asset supply','Location','best');

%saveas(gcf,'_figures/Figure2bReplication.png');