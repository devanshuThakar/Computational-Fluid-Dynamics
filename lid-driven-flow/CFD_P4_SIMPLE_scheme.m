%% CFD Project-4: To solve Lid-driven cavity flow.
%   Scheme: SIMPLE

clc;
clear variables;
close all;

method='SIMPLE';
tic

%% Defining the domain:
N = 100;                     % Number of control volumes in 'x' and 'y' direction.
dom_length = 1;
dh = dom_length/(N-1);
x = 0:dh:dom_length;            % X domain length
y = 0:dh:dom_length;            % Y domain length
nu = 0.01;                      % Kinematic viscocity.

% Under-relaxation factors:
URF = 0.8;
URF_p = 0.8;

%% Initializing the variables

%Final collocated variables
    u_final = zeros(N,N);
    v_final = zeros(N,N);
    
    %% Step-1: Guess the pressure field.
    p_final = ones(N,N);
    u_final(1,:) = 1;
    %%
    
% Staggered variables:
    u = zeros(N+1,N);
    u_stag = zeros(N+1,N);
    d_e = zeros(N+1,N);
    v = zeros(N,N+1);
    v_stag = zeros(N,N+1);
    d_n = zeros(N,N+1);
    p = ones(N+1,N+1);
    p_stag = ones(N+1,N+1);
    pc = zeros(N+1,N+1);
    b = zeros(N+1,N+1);
    u(1,:)= 2;

u_new = zeros(N+1,N);
v_new = zeros(N,N+1);
p_new = ones(N+1,N+1);
u_new(1,:)=1;

%% Step-2:  Solving the discritized FVM equations [Eq. 3 and 4]: 
err = 1;
iter = 0;
error_tol = 1e-6;   % Tolerance error.


while err > error_tol
    
    % X-mom equations - For Interior CVs
    for i = 2:N
        for j = 2:N - 1
            u_e = 0.5*(u(i,j) + u(i,j+1));
            u_w = 0.5*(u(i,j) + u(i,j-1));
            v_n = 0.5*(v(i-1,j) + v(i-1,j+1));
            v_s = 0.5*(v(i,j) + v(i,j+1));
            
            a_E = -0.5*u_e*dh + nu;
            a_W = 0.5*u_w*dh + nu;
            a_N = -0.5*v_n*dh + nu;
            a_S = 0.5*v_s*dh + nu;
            
            a_e = 0.5*u_e*dh - 0.5*u_w*dh + 0.5*v_n*dh - 0.5*v_s*dh + 4*nu;
            
            A_e = -dh;
            d_e(i,j) = A_e/a_e;
            
            u_stag(i,j) = (a_E*u(i,j+1) + a_W*u(i,j-1) + a_N*u(i-1,j) + a_S*u(i+1,j))/a_e + d_e(i,j)*(p(i,j+1) - p(i,j));
        end
    end
    
    % X-mom equations - For Boundary CVs
        u_stag(1,:) = 2 - u_stag(2,:);
        u_stag(N + 1,:) = -u_stag(N,:);
        u_stag(2:N,1) = 0;
        u_stag(2:N,N) = 0;
    
    % Y-mom equations - For Interior CVs
    for i = 2:N - 1
        for j = 2:N
            u_e = 0.5*(u(i,j) + u(i+1,j));
            u_w = 0.5*(u(i,j-1) + u(i+1,j-1));
            v_n = 0.5*(v(i-1,j) + v(i,j));
            v_s = 0.5*(v(i,j) + v(i+1,j));
            
            a_E = -0.5*u_e*dh + nu;
            a_W = 0.5*u_w*dh + nu;
            a_N = -0.5*v_n*dh + nu;
            a_S = 0.5*v_s*dh + nu;
            
            a_n = 0.5*u_e*dh - 0.5*u_w*dh + 0.5*v_n*dh - 0.5*v_s*dh + 4*nu;
            
            A_n = -dh;
            d_n(i,j) = A_n/a_n;
            
            v_stag(i,j) = (a_E*v(i,j+1) + a_W*v(i,j-1) + a_N*v(i-1,j) + a_S*v(i+1,j))/a_n + d_n(i,j)*(p(i,j) - p(i+1,j));
        end
    end
    
    % Y-mom equations - For Boundary CVs
        v_stag(:,1) = -v_stag(:,2);
        v_stag(:,N + 1) = -v_stag(:,N);
        v_stag(1,2:N) = 0;
        v_stag(N,2:N) = 0;
    
    % Again intializing the Pressure-correction term:
        pc(1:N+1,1:N+1)=0;
    
    %% Step-3: Solving for Continuity equation to get P' - For Interior CVs
    for i = 2:N
        for j = 2:N
            a_E = -d_e(i,j)*dh;
            a_W = -d_e(i,j-1)*dh;
            a_N = -d_n(i-1,j)*dh;
            a_S = -d_n(i,j)*dh;
            a_P = a_E + a_W + a_N + a_S;
            b(i,j) = -(u_stag(i,j) - u_stag(i,j-1))*dh + (v_stag(i,j) - v_stag(i-1,j))*dh;
            
            %% Finding pc using Gauss-Jacobi algorithm:
            pc(i,j) = (a_E*pc(i,j+1) + a_W*pc(i,j-1) + a_N*pc(i-1,j) + a_S*pc(i+1,j) + b(i,j))/a_P;
        end
    end
    
    %% Step-4: Correcting the Pressure field:
    for i = 2:N
        for j = 2:N
            p_new(i,j) = p(i,j) + URF_p*pc(i,j);
        end
    end
    
    % Continuity equation - Boundary CVs:
        % Implementing pressure boundary condtions: 
        %       do(p)/do(x) = 0 at west and east face. 
        %       do(p)/do(y)=0 at north and south faces.
        
        p_new(1,:) = p_new(2,:);
        p_new(N + 1,:) = p_new(N,:);
        p_new(:,1) = p_new(:,2);
        p_new(:,N + 1) = p_new(:,N);
    
    % Correcting the velocities:
    for i = 2:N
        for j = 2:N - 1
            u_new(i,j) = u_stag(i,j) + URF*d_e(i,j)*(pc(i,j+1) - pc(i,j));
        end
    end
    
    % X_mom equations - Boundary CVs.
        u_new(1,:) = 2 - u_new(2,:);
        u_new(N + 1,:) = -u_new(N,:);
        u_new(2:N,1) = 0;
        u_new(2:N,N) = 0;
    
    for i = 2:N - 1
        for j = 2:N
            v_new(i,j) = v_stag(i,j) + URF*d_n(i,j)*(pc(i,j) - pc(i+1,j));
        end
    end
    
    % Y_mom equations - Boundary CVs.
        v_new(:,1) = -v_new(:,2);
        v_new(:,N + 1) = -v_new(:,N);
        v_new(1,2:N) = 0;
        v_new(N,2:N) = 0;
            
    
    % Updating error value from residual term 'b':
    err = 0;
    for i = 2:N
        for j = 2:N
            err = err + abs(b(i,j));
        end
    end
    
    %% Step-5: Updating the values of variables and completing the iterations:
    u = u_new;
    v = v_new;
    p = p_new;
    iter = iter + 1;
    
    disp(['Iteration: ', num2str(iter),' Error : ', num2str(err)])
    
end
toc

% Mapping staggered variables to collocated variables:
for i = 1:N
    for j = 1:N
        u_final(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_final(i,j) = 0.5*(v(i,j) + v(i,j+1));
        % No change for pressure term as it's values are stored at collacated grid.
        p_final(i,j) = p(i,j);
    end
end


%% Plotting the results:
x_dom = ((1:N)-1).*dh;
y_dom = 1-((1:N)-1).*dh;
[X,Y] = meshgrid(x_dom,y_dom);

figure(1);
    contourf(X,Y,u_final, 21, 'LineStyle', 'none')
    colorbar
    colormap('jet')
    title(sprintf(' Contour plot of u in the domain for %sX%s',num2str(N),num2str(N)));
    xlabel('x')
    ylabel('y')
    saveas(gcf,sprintf('%s_Contour_plot_u_%sX%s_grid.png',method,num2str(N),num2str(N)));

figure(2);
    contourf(X,Y,v_final, 21, 'LineStyle', 'none')
    colorbar
    colormap('jet')
    title(sprintf(' Contour plot of v in the domain for %sX%s',num2str(N),num2str(N)));
    xlabel('x')
    ylabel('y')
    saveas(gcf,sprintf('%s_Contour_plot_v_%sX%s_grid.png',method,num2str(N),num2str(N)));

figure(3);
    quiver(X, Y, u_final, v_final, 3, 'r')
    xlabel('x')
    ylabel('y')
    axis equal
    grid on
    title(sprintf('Streamline plot for %sX%s',num2str(N),num2str(N)));
    saveas(gcf,sprintf('%s_Streamline_plot_%sX%s_grid.png',method,num2str(N),num2str(N)));

figure(4)
    plot(y(:),u_final(:,N/2));
    hold on;
    plot(y(:),v_final(:,N/2));
    xlim([0 1])
    title(sprintf('Plot of u,v at x=0.5 centerline for %sX%s',num2str(N),num2str(N)));
    legend('u','v');
    ylabel('u,v');
    xlabel('y');
    grid on;
    saveas(gcf,sprintf('%s_plot_uv_at_x_centerline_%sX%s_grid.png',method,num2str(N),num2str(N)));

figure(5);
    plot(x(:),u_final(N/2,:));
    hold on;
    plot(x(:),v_final(N/2,:));
    xlim([0 1])
    grid on;
    title(sprintf('Plot of u,v at y=0.5 centerline for %sX%s',num2str(N),num2str(N)));
    legend('u','v');
    ylabel('u,v');
    xlabel('x');
    saveas(gcf,sprintf('%s_plot_uv_at_y_centerline_%sX%s_grid.png',method,num2str(N),num2str(N)));

figure(6);
    plot(x(:),p_final(N/2,:));
    hold on;
    plot(y(:),p_final(:,N/2));
    grid on;
    title(sprintf('Plot of P at x=0.5,y=0.5 centerline for %sX%s',num2str(N),num2str(N)));
    xlim([0 1]);
    legend('at y=0.5','at x=0.5');
    ylabel('P');
    xlabel('x,y');
    saveas(gcf,sprintf('%s_plot_P_at_xy_centerline_%sX%s_grid.png',method,num2str(N),num2str(N)));
