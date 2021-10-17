clear all;
clc
%% Plot of residual with iteration for incresing grid points
load('residual_X_40_Y_40_ADI.mat')
grid_40 = residual;
load('residual_X_80_Y_80_ADI.mat');
grid_80 = residual;
load('residual_X_160_Y_160_ADI.mat')
grid_160 = residual;

loglog(grid_40(1:end), 'LineWidth',1)
hold on
loglog(grid_80(1:end), 'LineWidth',1)
hold on
loglog(grid_160(1:end),'LineWidth',1)
title('Plot of residual vs. iterations by ADI method');
legend('40 X 40', '80 X 80', '160 X 160');
ylabel('residual');
xlabel('Iteration');
grid on
saveas(gcf,'ADI_residual_Iteration_PartA_Q4.png')

%% Contour plot of phi for 160 X 160 grid by ADI
load('phi_num_160_160_ADI.mat')
figure(2)
contour(x,y,phi_numerical,50)
title('Contour plot of phi by ADI method');
saveas(gcf, 'Contour_by_ADI_for_160_PartA_Q4.png')

%% Plot of residual for row-sweep, column-sweep and ADI for 40 X 80 grid
load('residual_X_40_Y_80_row_sweep.mat')
row = residual;
load('residual_X_40_Y_80_column_sweep.mat');
column = residual;
load('residual_X_40_Y_80_ADI.mat')
ADI = residual;

figure(3)
loglog(row(1:end), 'LineWidth',1)
hold on
loglog(column(1:end), 'LineWidth',1)
hold on
loglog(ADI(1:end),'LineWidth',1)
title('Plot of residual vs. iterations');
legend('row sweep', 'column sweep', 'ADI');
ylabel('residual');
xlabel('Iteration');
grid on
saveas(gcf,'row_column_ADI_residual_Iteration_PartA_Q4.png')