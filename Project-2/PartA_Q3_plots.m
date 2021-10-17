%% Plot of residual for 80 X 80 grid
load('residual_X_80_Y_80_row_sweep.mat');
row_80 = residual; 

load('residual_X_80_Y_80_column_sweep.mat');
column_80 = residual;

p1 = loglog(row_80(1:end), 'LineWidth',1);
p1.Color(4) = 1;
hold on
p2=loglog(column_80(1:end), 'LineWidth',1);
p2.Color(4) = 0.75;
hold on
title(sprintf('Plot of residual vs. iterations'));
legend('row sweep', 'column sweep');
ylabel('residual');
xlabel('Iteration');
grid on
saveas(gcf, sprintf('X80_Y80_Residual_Iteration_PartA_Q3.png'))

%% Plot of residual for 40 X 80 grid
load('residual_X_40_Y_80_row_sweep.mat');
row_80 = residual; 

load('residual_X_40_Y_80_column_sweep.mat');
column_80 = residual;
figure(2);
p1 = loglog(row_80(1:end), 'LineWidth',1);
p1.Color(4) = 1;
hold on
p2=loglog(column_80(1:end), 'LineWidth',1);
p2.Color(4) = 0.75;
hold on
title(sprintf('Plot of residual vs. iterations'));
legend('row sweep', 'column sweep');
ylabel('residual');
xlabel('Iteration');
grid on
saveas(gcf, sprintf('X40_Y80_Residual_Iteration_PartA_Q3.png'))

%% Contour plots
load('phi_num_X80_Y80_row_sweep.mat')
x_=x;
y_=y;
phi_numerical_=phi_numerical;
figure(3)
contour(x_,y_,phi_numerical_,50)
title(sprintf('Contour plot of phi by row sweep for 80 X 80'));
saveas(gcf, sprintf('Row_sweep_Contour_80_PartA_Q3.png'))

load('phi_num_X40_Y80_row_sweep.mat')
figure(4)
contour(y,x,phi_numerical,50)
title(sprintf('Contour plot of phi by row sweep for 40 X 80'));
saveas(gcf, sprintf('Row_sweep_Contour_40X_80Y_PartA_Q3.png'))