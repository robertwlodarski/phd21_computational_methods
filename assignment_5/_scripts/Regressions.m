% Regression analysis 

% Prepare data & run regresison 
vLogK           = log(Results.vK);
vLogA           = log(Results.vA);
X               = [ones(length(vLogA),1) vLogA];

% Perform linear regression
b               = X \ vLogK;  
vLogKFitted     = X*b;
vEpsilon        = vLogK-vLogKFitted;

% Plot the figure
figure(1);
scatter(vLogA, vLogK,'filled','square','AlphaData',0.5); 
hold on; 
grid on;
plot(vLogA, vLogKFitted, 'r-', 'LineWidth', 2);
xlabel('log(A)');
ylabel('log(K)');
legend('Data', 'Fitted line');
hold off;
saveas(gcf,'_figures/FittedRegression.png');

% Plot the figure
figure(2);
plot(1:1:length(vLogK), vEpsilon,'Color','black','LineWidth', 1); 
hold on; 
grid on;
xlabel('T');
ylabel('$\epsilon$','Interpreter','latex','FontSize',14);
hold off;
saveas(gcf,'_figures/ErrorTime.png');

% Plot the figure
figure(3);
scatter(vLogA, vEpsilon,'filled','square','AlphaData',0.5); 
grid on;
xlabel('log(A)');
ylabel('$\epsilon$','Interpreter','latex','FontSize',14);
saveas(gcf,'_figures/ErrorA.png');

% Plot the figure
figure(4);
plot(1:1:length(vLogK), vLogK,'b-', 'LineWidth', 1); 
hold on; 
grid on;
plot(1:1:length(vLogK), vLogKFitted, 'r-', 'LineWidth', 1);
xlabel('T');
ylabel('log(K)');
legend('Data', 'Fitted line');
hold off;
saveas(gcf,'_figures/FittedK.png');

% Plot the figure
figure(5);
plot(1:1:500, vLogK(1:1:500),'b-', 'LineWidth', 1.5); 
hold on; 
grid on;
plot(1:1:500, vLogKFitted(1:1:500), 'r-', 'LineWidth', 1.5);
xlabel('T');
ylabel('log(K)');
legend('Data', 'Fitted line');
hold off;
saveas(gcf,'_figures/FittedK500.png');