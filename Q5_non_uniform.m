clear all;
clc;
% Defining Parameters
n=11;
L=1;
rho=1;
u=1;
phi_o=0;
phi_L=1;
Pe=50;
Gamma = rho*u*L/Pe;

%Since the given n also includes the
X_exact = linspace(0,1,100);
X = zeros(1,n);
r = 0.7; % expansion factor
delta_x_max = L*(r-1)/(power(r,n)-1);  % For non-unifrom contracting grid
%Assigning X
X(1)=0;
for i=2:n
    X(i)=X(i-1) + delta_x_max*power(r,i-2);
end
n_interi=n-2;
X_interi = X(2:n_interi);
Ap=zeros(1,n_interi);
Ae=zeros(1,n_interi);
Aw=zeros(1,n_interi);
for i=2:n-1
    j=i-1;
    %% Comment below three lines to change second derivative from CDS to FDS
%     Diffusion='FDS';
%     Ap(j) = rho*u/(X(i)-X(i-1)) + (Gamma*(X(i+1)-X(i-1)))/((X(i+1)-X(i))^2 * (X(i)-X(i-1)));
%     Ae(j) = -Gamma*(X(i)-X(i-1))/( (X(i+1)-X(i))^2 * (X(i)-X(i-1)) );
%     Aw(j) = -rho*u/(X(i)-X(i-1)) - (Gamma*(X(i+1)-X(i)))/((X(i+1)-X(i))^2 * (X(i)-X(i-1)));
    %% Comment above three lines to cahnge from FDS to CDS
    Diffusion='CDS';
    Ap(j) = rho*u/(X(i)-X(i-1)) + (Gamma*(X(i+1)-X(i-1)))/(0.5*(X(i+1)-X(i-1))*(X(i+1)-X(i)) * (X(i)-X(i-1)));
    Ae(j) = -Gamma*(X(i)-X(i-1))/(0.5*(X(i+1)-X(i-1))*(X(i+1)-X(i)) * (X(i)-X(i-1)));
    Aw(j) = -rho*u/(X(i)-X(i-1)) - (Gamma*(X(i+1)-X(i)))/(0.5*(X(i+1)-X(i-1))*(X(i+1)-X(i)) * (X(i)-X(i-1)));
end
B = zeros(1,n_interi);
B(1)=-Aw(1)*phi_o;
B(n_interi)=-Ae(n_interi)*phi_L;

y_exact=phi_analytical(Pe,L,phi_o,phi_L,X_exact);
y_num = zeros(n,1);
y_num(1)=phi_o;
y_num(n)=phi_L;
y_num(2:n-1)=ThomasAlgo(Ap,Aw(2:n_interi),Ae(1:n_interi-1),B);

f1 = figure();
plot(X_exact,y_exact);
hold on;
plot(X,y_num);
title('Plot of \phi vs x');
legend('Exact', 'Numerical');
xticks(0:0.1:1);
yticks(0:0.1:1);
set(gca,'XLim',[0 1], 'YLim', [0 1]);
xlabel('x');
ylabel('\phi');
grid on;
saveas(gcf, sprintf('NonUniform_%iNodes_BDS_in_Covection_%s_in_diffusion_%i_Pe_phi_vs_x.png', n, Diffusion, Pe))