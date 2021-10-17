clear all;
clc
method='GaussSeidel';

load(sprintf('residual_40_%s.mat',method))
grid_40 = residual;
load(sprintf('residual_80_%s.mat',method));
grid_80 = residual;
load(sprintf('residual_160_%s.mat',method))
grid_160 = residual;

loglog(grid_40(1:end), 'LineWidth',1)
hold on
loglog(grid_80(1:end), 'LineWidth',1)
hold on
loglog(grid_160(1:end),'LineWidth',1)
title(sprintf('Plot of residual vs. iterations for %s method',method));
legend('40 X 40', '80 X 80', '160 X 160');
ylabel('residual');
xlabel('Iteration');
grid on
saveas(gcf, sprintf('%s_residual_Iteration_PartA_Q2.png',method))

load(sprintf('phi_num_160_%s.mat',method))
figure(2)
contour(x,y,phi_numerical,50)
title(sprintf('Contour plot of phi by %s method',method));
saveas(gcf, sprintf('%s_Contour_160_PartA_Q2.png',method))

figure(3)
contourf(x,y,phi_numerical,50)
title(sprintf('Contour plot of phi by %s method',method));
saveas(gcf, sprintf('%s_Contourf_160_PartA_Q2.png', method))