%% Plotting distributions

% Find maximum index (say a

figure(11);
bar(Grids.vGridA2,Results.mMarginalWealth(1,:),'blue')
hold on;
bar(Grids.vGridA2,Results.mMarginalWealth(50,:),'red')