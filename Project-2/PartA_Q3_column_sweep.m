clear all;
clc

N=80; %Number of grid points in x
M =80; %Number of grid points in y 
phi_num(:,:,2) = zeros(N,M); phi_num(:,:,1) = zeros(N,M);
x = linspace(0,1,N);
y = linspace(0,1,M);
delta_x = 1/(N-1); delta_y=1/(M-1);

phi_num(1,:,1) = phi(x(1)*ones(1,M), y); phi_num(1,:,2) = phi_num(1,:,1);
phi_num(N,:,1) = phi(x(N)*ones(1,M), y); phi_num(N,:,2) = phi_num(N,:,1);
phi_num(:,1,1) = phi(x,y(1)*ones(1,N)); phi_num(:,1,2) = phi_num(:,1,1);
phi_num(:,M,1) = phi(x,y(M)*ones(1,N)); phi_num(:,M,2) = phi_num(:,M,1);

source_arr = zeros(N,M);
for i=2:N-1
    for j=2:M-1
        source_arr(i,j) = source(x(i), y(j));
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
    % column sweep, index i runs from row 2 to row M-1
    for j=2:M-1
        Ap = zeros(N-2,1)-2*(1+ (delta_x^2)/(delta_y^2));
        Ae = zeros(N-2-1,1)+1;
        Aw = zeros(N-2-1,1)+1;
        B = source_arr(2:N-1,j)*delta_x^2 - (phi_num(2:N-1,j+1,prev) + phi_num(2:N-1,j-1,prev))*(delta_x^2)/(delta_y^2);
        
        % In B -need to consider the boundary values for first and last cell
        B(1) = B(1) - phi_num(1,j,curr); % Boundary can either be current or prev
        B(N-2) = B(N-2)-phi_num(N,j,curr);
        phi_num(2:N-1,j,curr) = transpose(ThomasAlgo(Ap,Aw,Ae,B));
        %phi_num(2:N-1,j,prev)=phi_num(2:N-1,j,curr); 
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
save(sprintf('phi_num_X%i_Y%i_column_sweep',N,M), 'phi_numerical', 'x', 'y');
save(sprintf('residual_X_%i_Y_%i_column_sweep', N,M), 'residual');