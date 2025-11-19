% Transitional dynamics plots 

%% A behaviours

figure(1);
plot(1:1:50,Paths.v1(1:50),'LineWidth',2.5);
grid on;
hold on;
plot(1:1:50,Paths.v2(1:50),'LineWidth',2.5);
plot(1:1:50,Paths.v3(1:50),'LineWidth',2.5);
legend('Permanent positive: A = 1.10 $\forall \; t \geq 1$',...
       'One-off negative: A = 1 - 0.15 \times \mathbf{1}($t=1)$',...
       'Transitory negative: $\log A_{t+1}=0.9^t \times \log 0.85 \; \forall t$','interpreter','latex','FontSize',10, 'Location','best');
saveas(gcf,'_figures/ShockTypes.png');

%% K behaviours

figure(2);
plot(1:1:Parameters.pT,Results.vK,'LineWidth',2.5);
grid on;
hold on;
plot(1:1:Parameters.pT,Results2.vK,'LineWidth',2.5);
plot(1:1:Parameters.pT,Results3.vK,'LineWidth',2.5);
plot(1:1:Parameters.pT,repmat(Results.vK0,Parameters.pT,1),'--','LineWidth',2.5);
plot(1:1:Parameters.pT,repmat(Results.vKT,Parameters.pT,1),'--','LineWidth',2.5);
xlabel('$t$','interpreter','latex','FontSize',14);
ylabel('$K$','interpreter','latex','FontSize',14);
legend('Permanent positive: A = 1.10 $\forall \; t \geq 1$',...
       'One-off negative: A = 1 - 0.15 \times \mathbf{1}($t=1)$',...
       'Transitory negative: $\log A_{t+1}=0.9^t \times \log 0.85 \; \forall t$', '$A=1$ SS', '$A=1.1$ SS','interpreter','latex','FontSize',10, 'Location','best');
saveas(gcf,'_figures/CapitalTransitions.png');

%% K behaviour: PE vs GE

figure(3);
plot(1:1:Parameters.pT,Results.vK,'LineWidth',2.5);
grid on;
hold on;
plot(1:1:Parameters.pT,ResultsPE.vK,'LineWidth',2.5);
plot(1:1:Parameters.pT,repmat(Results.vK0,Parameters.pT,1),'--','LineWidth',2.5);
plot(1:1:Parameters.pT,repmat(Results.vKT,Parameters.pT,1),'--','LineWidth',2.5);
xlabel('$t$','interpreter','latex','FontSize',14);
ylabel('$K$','interpreter','latex','FontSize',14);
legend('General equilbrium','Partial equilibrium','$A=1$ SS', '$A=1.1$ SS','interpreter','latex','FontSize',10, 'Location','best')
saveas(gcf,'_figures/CapitalTransitionsEq.png');

%% Gini

figure(4);
plot(1:1:Parameters.pT,Results.vGini,'LineWidth',2.5);
grid on;
hold on;
plot(1:1:Parameters.pT,Results2.vGini,'LineWidth',2.5);
plot(1:1:Parameters.pT,Results3.vGini,'LineWidth',2.5);
xlabel('$t$','interpreter','latex','FontSize',14);
ylabel('Gini','interpreter','latex','FontSize',14);
legend('Permanent positive: A = 1.10 $\forall \; t \geq 1$',...
       'One-off negative: A = 1 - 0.15 \times \mathbf{1}($t=1)$',...
       'Transitory negative: $\log A_{t+1}=0.9^t \times \log 0.85 \; \forall t$','interpreter','latex','FontSize',10, 'Location','best');
saveas(gcf,'_figures/GiniTransitions.png');

%% Gini behaviour: PE vs GE

figure(5);
plot(1:1:Parameters.pT,Results.vGini,'LineWidth',2.5);
grid on;
hold on;
plot(1:1:Parameters.pT,ResultsPE.vGini,'LineWidth',2.5);
xlabel('$t$','interpreter','latex','FontSize',14);
ylabel('Gini','interpreter','latex','FontSize',14);
legend('General equilbrium','Partial equilibrium','interpreter','latex','FontSize',10, 'Location','best')
saveas(gcf,'_figures/GiniTransitionsEq.png');