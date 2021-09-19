clear all;
clc;
% Defining Parameters
n=41;
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
Ap = ones(1,n_interior)*(4*Gamma);
Aw = ones(1,n_interior-1)*-(rho*u*delta_x + 2*Gamma);
Ae = ones(1,n_interior-1)*(rho*u*delta_x - 2*Gamma);

B = zeros(1,n_interior);
B(1) = (rho*u*delta_x + 2*Gamma)*phi_o;
B(n_interior) = -(rho*u*delta_x - 2*Gamma)*phi_L;

y_exact=phi_analytical(Pe,L,phi_o,phi_L,X_exact);
y_num = zeros(n,1);
y_num(1)=phi_o;
y_num(n)=phi_L;
y_num(2:n-1)=ThomasAlgo(Ap,Aw,Ae,B);

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
saveas(gcf, sprintf('Uniform_Nodes_CDSinCovection_%i_Pe_%i_phi_vs_x.png', n, Pe))

function y=ThomasAlgo(Ap,Aw,Ae,B)
n=length(Ap);
y=zeros(n,1);
%Decomposition
for i=1:n-1
    ratio = Aw(i)/Ap(i);
    Ap(i+1)=Ap(i+1) - Ae(i)*ratio;
    Aw(i) = ratio;
end

%Forward Substitution
for i=2:n
    B(i)=B(i)- B(i-1)*Aw(i-1);
end

%Backward Substitution
for i=n:-1:1
if(i==n)
    y(i)=B(i)/Ap(i);
    continue;
end
y(i) = (B(i) - y(i+1)*Ae(i))/Ap(i);
end
end

function y=phi_analytical(Pe,L,phi_o,phi_L, X)
y = (exp(X*Pe/L) - 1)/(exp(Pe) - 1);
y = y*(phi_L - phi_o);
y=y + phi_o;
end
