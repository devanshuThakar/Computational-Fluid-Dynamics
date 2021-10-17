clear all;
clc

N=40; %Number of grid points
phi_num(:,:,2) = zeros(N); phi_num(:,:,1) = zeros(N);
x = linspace(0,1,N);
delta_x = 1/(N-1);
y = linspace(0,1,N);

phi_num(1,:,1) = phi(x(1)*ones(1,length(x)), y); phi_num(1,:,2) = phi_num(1,:,1);
phi_num(N,:,1) = phi(x(N)*ones(1,length(x)), y); phi_num(N,:,2) = phi_num(N,:,1);
phi_num(:,1,1) = phi(x,y(1)*ones(1,length(x))); phi_num(:,1,2) = phi_num(:,1,1);
phi_num(:,N,1) = phi(x,y(N)*ones(1,length(x))); phi_num(:,N,2) = phi_num(:,N,1);

source_arr = zeros(N);
for i=2:N-1
    for j=2:N-1
        source_arr(i,j) = source(x(i), y(j));
    end
end

done = false;
iter=0;
tolerance = 1e-05;
residual = [];

while ~done && iter < 1000000
    not_converged = [];
%     done = true;
    if mod(iter,2)==0
        curr=2;
        prev=1;
    else
        curr=1;
        prev=2;
    end
    resi = 0;
    for i=2:N-1
        for j=2:N-1
            phi_num(i,j,curr) = (phi_num(i+1,j,prev) + phi_num(i-1,j,prev) + phi_num(i,j+1,prev)+phi_num(i,j-1,prev) - source_arr(i,j)*delta_x^2)/4;
        end
    end
    
    resi = source_arr(2:N-1,2:N-1)*delta_x^2 - (phi_num(3:N,2:N-1,curr) +  phi_num(1:N-2,2:N-1,curr) + phi_num(2:N-1,3:N,curr) + phi_num(2:N-1,1:N-2,curr) - 4*phi_num(2:N-1,2:N-1,curr));
    resi = norm(resi);
    if(resi < tolerance)
        done=true;
    end
    iter=iter+1;
    residual(iter) = resi;
end

phi_numerical = phi_num(:,:,curr);
% save(sprintf('residual_%i_Jacobi', N), 'residual');
% save(sprintf('phi_num_%i_Jacobi',N), 'phi_numerical', 'x', 'y');