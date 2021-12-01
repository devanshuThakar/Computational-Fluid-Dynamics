% load(sprintf('%s_%sX%s.mat',method,num2str(N-2),num2str(N-2)));

% plot(y(:),u_curr(:,N/2+1),'-o','MarkerSize',5);
figure(1)
plot(y(:),u_curr(:,N/2+1));
hold on;
plot(y(:),v_curr(:,N/2+1));
xlim([0 1])
title(sprintf('Plot of u,v at x=0.5 centerline for %sX%s',num2str(N-2),num2str(N-2)));
legend('u','v');
ylabel('u,v');
xlabel('y');
grid on;
saveas(gcf,sprintf('%s_plot_uv_at_x_centerline_%sX%s_grid.png',method,num2str(N-2),num2str(N-2)));

figure(2);
plot(x(:),u_curr(N/2+1,:));
hold on;
plot(x(:),v_curr(N/2+1,:));
xlim([0 1])
grid on;
title(sprintf('Plot of u,v at y=0.5 centerline for %sX%s',num2str(N-2),num2str(N-2)));
legend('u','v');
ylabel('u,v');
xlabel('x');
saveas(gcf,sprintf('%s_plot_uv_at_y_centerline_%sX%s_grid.png',method,num2str(N-2),num2str(N-2)));

figure(3);
plot(x(:),P_curr(N/2+1,:));
hold on;
plot(y(:),P_curr(:,N/2+1));
grid on;
title(sprintf('Plot of P at x=0.5,y=0.5 centerline for %sX%s',num2str(N-2),num2str(N-2)));
xlim([0 1]);
legend('at y=0.5','at x=0.5');
ylabel('P');
xlabel('x,y');
saveas(gcf,sprintf('%s_plot_P_at_xy_centerline_%sX%s_grid.png',method,num2str(N-2),num2str(N-2)));
