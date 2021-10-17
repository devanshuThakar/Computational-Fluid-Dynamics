% Implementation of Crank Nicolson scheme with ADI to solve time dependent 2D problems %
clear all;
clc

N=40; %Number of grid points in x
M =40; %Number of grid points in y 
ts=0;
te=10;
delta_t = 0.01;
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


done = false;
tolerance = 1e-05;
error_arr = [];
residual = [];
iter=0;

phi_temp = zeros(N,M);

while ~done && t < te
    if mod(iter,2)==0
        curr=2;
        prev=1;
    else
        curr=1;
        prev=2;
    end
    %phi_temp = zeros(N,M);
     %Do a row sweep
     % row sweep, index i runs from row 2 to row N-1
     for i=2:N-1
         Ap = zeros(M-2,1)-(4+(4*delta_x^2/delta_t));
         Ae = zeros(M-2-1,1)+1;
         Aw = zeros(M-2-1,1)+1;
         B = 2*source_arr(i,2:M-1)*delta_x^2 ...
             - (4*delta_x^2/delta_t)*phi_num(i,2:M-1,prev) ... 
             - 2*phi_num(i+1,2:M-1,prev) - 2*phi_num(i-1,2:M-1,prev) ...
             - phi_num(i,1:M-2,prev) - phi_num(i,3:M,prev) + 4*phi_num(i,2:M-1,prev);

         % In B -need to consider the boundary values for first and last cell
         B(1) = B(1) - phi_num(i,1,curr); % Boundary can either be current or prev
         B(M-2) = B(M-2)-phi_num(i,M,curr);
         phi_temp(i,2:M-1) = transpose(ThomasAlgo(Ap,Aw,Ae,B));
     end
     
     % Do a column sweep
     % column sweep, index i runs from row 2 to row M-1
     for j=2:M-1
         Ap = zeros(N-2,1)-(4+(4*delta_x^2/delta_t));
         Ae = zeros(N-2-1,1)+1;
         Aw = zeros(N-2-1,1)+1;
         B = 2*source_arr(2:N-1,j)*delta_x^2  - (4*(delta_x^2)/delta_t)*phi_temp(2:N-1,j)...
             - phi_num(3:N,j,prev) - phi_num(1:N-2,j,prev) - phi_num(1:N-2,j+1,prev) - phi_num(1:N-2,j-1,prev) + 4*phi_num(1:N-2,j,prev)...
             - phi_temp(2:N-1,j+1)-phi_temp(2:N-1,j-1);

         % In B -need to consider the boundary values for first and last cell
         B(1) = B(1) - phi_num(1,j,curr); % Boundary can either be current or prev
         B(N-2) = B(N-2)-phi_num(N,j,curr);
         phi_num(2:N-1,j,curr) = transpose(ThomasAlgo(Ap,Aw,Ae,B));
      end
    
    if(iter==0)
        phi_num(1,:,1) = phi(x(1)*ones(1,length(x)), y); phi_num(1,:,2) = phi_num(1,:,1); phi_temp(1,:)=phi_num(1,:,1);
        phi_num(N,:,1) = phi(x(N)*ones(1,length(x)), y); phi_num(N,:,2) = phi_num(N,:,1); phi_temp(N,:)=phi_num(N,:,1); 
        phi_num(:,1,1) = phi(x,y(1)*ones(1,length(y))); phi_num(:,1,2) = phi_num(:,1,1); phi_temp(:,1)=phi_num(:,1,1);
        phi_num(:,M,1) = phi(x,y(M)*ones(1,length(y))); phi_num(:,M,2) = phi_num(:,M,1); phi_temp(:,M)=phi_num(:,M,1);
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