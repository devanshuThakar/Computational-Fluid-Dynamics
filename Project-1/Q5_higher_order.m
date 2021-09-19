clear all;
clc;
% Defining Parameters
n=11;
L=1;
rho=1;
u=1;
phi_o=0;
phi_L=1;
Pe=18;
Gamma = rho*u*L/Pe;

%Since the given n also includes the
X_exact = linspace(0,1,100);
X = linspace(0,1,n);
delta_x = (1-0)/(n-1);  % For uniform spacing
% Numerically we need to be solve only for n-2 interior nodes
n_interior=n-2;
A = zeros(n_interior, n_interior);
A(1,1:2) = [-(rho*u*delta_x + 2*Gamma), Gamma];
A(2,1:4) = [-(12*rho*u*delta_x + 16*Gamma),12*rho*u*delta_x + 30*Gamma, -16*Gamma, Gamma];
A(n_interior-1,n_interior-3:n_interior) = [Gamma, -(12*rho*u*delta_x + 16*Gamma),12*rho*u*delta_x + 30*Gamma, -16*Gamma];
A(n_interior, n_interior-1:n_interior) = [Gamma+rho*u*delta_x,-(rho*u*delta_x) + 2*Gamma];

for j=3:n_interior-2
    A(j,1+j-3)=Gamma;
    A(j,2+j-3)= -(12*rho*u*delta_x + 16*Gamma);
    A(j,3+j-3)= 12*rho*u*delta_x + 30*Gamma;
    A(j,4+j-3) = -16*Gamma;
    A(j,5+j-3) = Gamma;
end

B = zeros(1,n_interior);
B(1) = -(Gamma + rho*u*delta_x)*phi_o;
B(2) = -Gamma*phi_o;
B(n_interior-1) = -Gamma*phi_L;
B(n_interior) = -(Gamma*phi_L);

y_exact=phi_analytical(Pe,L,phi_o,phi_L,X_exact);
y_num = zeros(n,1);
y_num(1)=phi_o;
y_num(n)=phi_L;
y_num(2:n-1)=GaussEli(A,B);

f1 = figure();
plot(X_exact,y_exact);
hold on;
plot(X,y_num);
title('Plot of \phi vs x');
legend('Exact', 'Numerical');
xticks(0:0.1:1);
yticks(0:0.1:1);
% set(gca,'XLim',[0 1], 'YLim', [0 1]);
xlabel('x');
ylabel('\phi');
grid on;
saveas(gcf, sprintf('Higher_order_Nodes_%i_Pe_%i_phi_vs_x.png', n, Pe))