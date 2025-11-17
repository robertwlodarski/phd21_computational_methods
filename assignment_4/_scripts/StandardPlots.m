% Transitional dynamics plots 

%% K behaviour

figure(1);
plot(1:1:Parameters.pT,Results.vK,'LineWidth',2.5);
grid on;
xlabel('$t$','interpreter','latex','FontSize',14);
ylabel('$K$','interpreter','latex','FontSize',14);
saveas(gcf,'_figures/CapitalTransitions.png');