clear all;
clc

N=40; %Number of grid points in x
M =40; %Number of grid points in y 
ts=0;
te=10;
delta_t = 0.0001;
phi_num(:,:,2) = zeros(N,M); phi_num(:,:,1) = zeros(N,M);
t=ts;
x = linspace(0,1,N);
y = linspace(0,1,M);
delta_x = 1/(N-1); delta_y=1/(M-1);

phi_ss = zeros(N,M);
for i=1:N
    for j=1:M
        phi_ss(i,j) = phi(x(i),y(j));
    end
end
source_arr = zeros(N,M);
for i=2:N-1
    for j=2:M-1
        source_arr(i,j) = source(x(i), y(j));
    end
end

phi_num(1,:,1) = phi(x(1)*ones(1,length(x)), y); phi_num(1,:,2) = phi_num(1,:,1);
phi_num(N,:,1) = phi(x(N)*ones(1,length(x)), y); phi_num(N,:,2) = phi_num(N,:,1);
phi_num(:,1,1) = phi(x,y(1)*ones(1,length(y))); phi_num(:,1,2) = phi_num(:,1,1);
phi_num(:,M,1) = phi(x,y(M)*ones(1,length(y))); phi_num(:,M,2) = phi_num(:,M,1);
iter=1;
t=t+delta_t;

done = false;
tolerance = 1e-05;
error_arr = [];
residual = [];

while ~done && t < te
    if mod(iter,2)==0
        curr=2;
        prev=1;
    else
        curr=1;
        prev=2;
    end
    
    
    for i=2:N-1
        for j=2:M-1
            phi_num(i,j,curr) = phi_num(i,j,prev) + ...
                                (phi_num(i+1,j,prev) + phi_num(i-1,j,prev) - 2*phi_num(i,j,prev))*delta_t/(delta_x^2) + ...
                                (phi_num(i,j+1,prev) + phi_num(i,j-1,prev) - 2*phi_num(i,j,prev))*delta_t/(delta_y^2) - source_arr(i,j)*delta_t;
        end
    end
    resi = phi_num(:,:,curr) - phi_ss;
    resi = norm(resi);
    error = phi_num(:,:,curr)-phi_num(:,:,prev);
    error = norm(error);
    iter=iter+1;
    residual(iter) = resi;
    error_arr(iter)=error;
    t=t+delta_t;
    
    if(resi<tolerance)
        done=true;
    end
end