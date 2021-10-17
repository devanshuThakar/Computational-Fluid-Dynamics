clear all;
clc; 

N=80; %Number of grid points
phi_num = zeros(N);
x = linspace(0,1,N);
delta_x = 1/(N-1);
y = linspace(0,1,N);

% Note row index i represents y; column index j represents x
phi_num(1,:) = phi(x,y(1)*ones(1,length(x)));
phi_num(N,:) = phi(x,y(N)*ones(1,length(x)));
phi_num(:,1) = phi(x(1)*ones(1,length(x)), y);
phi_num(:,N) = phi(x(N)*ones(1,length(x)), y);

N_interi = N-2; % Number of interior grid points

A = zeros(N_interi*N_interi);
B = zeros(N_interi*N_interi, 1);

source_arr = zeros(N,N);
phi_anal = phi_num;
for i=2:N-1
    for j=2:N-1
        source_arr(i,j) = source(x(i), y(j));
        phi_anal(i,j) = phi(x(i),y(j));
    end
end

%Excluding the boundaries for A
for i=3:N-2
    for j=3:N-2
        % current row in A
        curr_row = (i-2)*N_interi+(j-1);
        B(curr_row) = source_arr(i,j)*delta_x^2;
        A((i-1-2)*N_interi+(j-1), (i-2)*N_interi+(j-1))=1;
        A((i+1-2)*N_interi+(j-1), (i-2)*N_interi+(j-1))=1;
        A((i-2)*N_interi+(j-1), (i-2)*N_interi+(j-1-1))=1;
        A((i-2)*N_interi+(j-1),(i-2)*N_interi+(j+1-1))=1;
        A((i-2)*N_interi+(j-1), (i-2)*N_interi+(j-1))=-4;
    end
end

%For boundary rows
for j=3:N-2
    i=2;
    curr_row = (i-2)*N_interi+(j-1);
    right = source_arr(i,j)*delta_x^2 - phi_num(1,j);
    A(curr_row,(i-2)*N_interi+(j-1-1))=1;
    A(curr_row,(i-2)*N_interi+(j+1-1))=1;
    A((i+1-2)*N_interi+(j-1),curr_row)=1;
    A(curr_row,curr_row)=-4;
    B(curr_row) = right;
    
    i=N-1;
    curr_row = (i-2)*N_interi+(j-1);
    right = source_arr(i,j)*delta_x^2 - phi_num(N,j);
    A(curr_row,(i-2)*N_interi+(j-1-1))=1;
    A(curr_row,(i-2)*N_interi+(j+1-1))=1;
    A((i+1-2)*N_interi+(j-1),curr_row)=1;
    A(curr_row,curr_row)=-4;
    B(curr_row) = right;
end

%For boundary columns
for i=3:N-2
    j=2;
    curr_row = (i-2)*N_interi+(j-1);
    right = source(x(i),y(i))*delta_x^2 - phi_num(i,1);
    A(curr_row,(i-2)*N_interi+(j+1-1))=1;
    A((i+1-2)*N_interi+(j-1),curr_row)=1;
    A((i-1-2)*N_interi+(j-1),curr_row)=1;
    A(curr_row,curr_row)=-4;
    B(curr_row) = right;
    
    j=N-1;
    curr_row = (i-2)*N_interi+(j-1);
    right = source(x(i),y(i))*delta_x^2 - phi_num(j,N);
    A(curr_row,(i-2)*N_interi+(j-1-1))=1;
    A((i+1-2)*N_interi+(j-1),curr_row)=1;
    A((i-1-2)*N_interi+(j-1),curr_row)=1;
    A(curr_row,curr_row)=-4;
    B(curr_row) = right;
end

%For corner points
curr_row = (2-2)*N_interi+(2-1);
i=2;j=2;
A(curr_row,translate(i,j+1,N))=1; A(translate(i+1,j,N),curr_row)=1; A(curr_row,curr_row)=-4;
B(curr_row) = source_arr(i,j)*delta_x^2 - phi_num(1,2) - phi_num(2,1);

curr_row = (N-1-2)*N_interi+(2-1);
i=N-1;j=2;
A(curr_row,translate(i,j+1,N))=1; A(translate(i-1,j,N),curr_row)=1; A(curr_row,curr_row)=-4;
B(curr_row) = source_arr(i,j)*delta_x^2 - phi_num(N,2) - phi_num(N-1,1);

curr_row = (2-2)*N_interi+(N-1-1);
i=2;j=N-1;
A(curr_row,translate(i,j-1,N))=1; A(translate(i+1,j,N),curr_row)=1; A(curr_row,curr_row)=-4;
B(curr_row) = source_arr(i,j)*delta_x^2 - phi_num(1,N-1) - phi_num(2,N);

curr_row = (N-1-2)*N_interi+(N-1-1);
i=N-1;j=N-1;
A(curr_row,translate(i,j-1,N))=1; A(translate(i-1,j,N),curr_row)=1; A(curr_row,curr_row)=-4;
B(curr_row) = source_arr(i,j)*delta_x^2 - phi_num(N,N-1) - phi_num(N-1,N);

t_start = cputime;
phi_GaussEli = GaussEli(A,B);
t_end = cputime -t_start;

for i=2:N-1
    for j=2:N-1
        phi_num(i,j) = phi_GaussEli(translate(i,j,N));
    end
end

resi = source_arr(2:N-1,2:N-1)*delta_x^2 - (phi_num(3:N,2:N-1) +  phi_num(1:N-2,2:N-1) + phi_num(2:N-1,3:N) + phi_num(2:N-1,1:N-2) - 4*phi_num(2:N-1,2:N-1));
resi = norm(resi);

function f=translate(i,j,N)
    N_interi=N-2;
    f=(i-2)*N_interi+(j-1);
end
