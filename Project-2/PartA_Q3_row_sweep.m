clear all;
clc

N=40; %Number of grid points in x
M =160; %Number of grid points in y 
phi_num(:,:,2) = zeros(N,M); phi_num(:,:,1) = zeros(N,M);
x = linspace(0,1,N);
y = linspace(0,1,M);
delta_x = 1/(N-1); delta_y=1/(M-1);

phi_num(1,:,1) = phi(x(1)*ones(1,M), y); phi_num(1,:,2) = phi_num(1,:,1);
phi_num(N,:,1) = phi(x(N)*ones(1,M), y); phi_num(N,:,2) = phi_num(N,:,1);
phi_num(:,1,1) = phi(x,y(1)*ones(1,N)); phi_num(:,1,2) = phi_num(:,1,1);
phi_num(:,M,1) = phi(x,y(M)*ones(1,N)); phi_num(:,M,2) = phi_num(:,M,1);

phi_anal = phi_num(:,:,1);
source_arr = zeros(N,M);
for i=2:N-1
    for j=2:M-1
        source_arr(i,j) = source(x(i), y(j));
        phi_anal(i,j) = phi(x(i),y(j));
    end
end

done = false;
iter=0;
tolerance = 1e-05;
residual = [];

while ~done && iter < 100000
    if mod(iter,2)==0
        curr=2;
        prev=1;
    else
        curr=1;
        prev=2;
    end
    % row sweep, index i runs from row 2 to row N-1
    for i=2:N-1
        Ap = zeros(M-2,1)-2*(1+(delta_y^2/delta_x^2));
        Ae = zeros(M-2-1,1)+1;
        Aw = zeros(M-2-1,1)+1;
        B = source_arr(i,2:M-1)*(delta_y^2) - (phi_num(i+1,2:M-1,prev) + phi_num(i-1,2:M-1,prev))*delta_y^2/delta_x^2;
        
        % In B -need to consider the boundary values for first and last cell
        B(1) = B(1) - phi_num(i,1,curr); % Boundary can either be current or prev
        B(M-2) = B(M-2)-phi_num(i,M,curr);
        phi_num(i,2:M-1,curr) = transpose(ThomasAlgo(Ap,Aw,Ae,B));
        %phi_num(i,2:M-1,prev)=phi_num(i,2:M-1,curr); 
    end
    
    resi = source_arr(2:N-1,2:M-1) - ... 
        (phi_num(3:N,2:M-1,curr) +  phi_num(1:N-2,2:M-1,curr)- 2*phi_num(2:N-1,2:M-1,curr))/delta_x^2 ...
        - (phi_num(2:N-1,3:M,curr) + phi_num(2:N-1,1:M-2,curr) - 2*phi_num(2:N-1,2:M-1,curr))/delta_y^2;
    
    resi = norm(resi);
    if(resi < tolerance)
        done=true;
    end
    iter=iter+1;
    residual(iter) = resi;
end

phi_numerical = phi_num(:,:,curr);
% save(sprintf('phi_num_X%i_Y%i_row_sweep',N,M), 'phi_numerical', 'x', 'y');
save(sprintf('residual_X_%i_Y_%i_row_sweep', N,M), 'residual');