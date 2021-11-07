load('error_UDS.mat')
load('error_QUICK.mat')
load('error_CDS.mat')

p1=loglog(error_CDS, 'LineWidth',1);
p1.Color(4) = 1;
hold on
p2=loglog(error_UDS, 'LineWidth',1);
p2.Color(4) = 0.75;
hold on
loglog(error_QUICK,'LineWidth',1)
title(sprintf('Plot of error vs. iterations for Gamma=0.01 and 80X80 grid'));
legend('CDS', 'UDS', 'QUICK');
ylabel('error');
xlabel('Iteration');
grid on
saveas(gcf, sprintf('error_vs_iteration_0.01_80X80.png'));