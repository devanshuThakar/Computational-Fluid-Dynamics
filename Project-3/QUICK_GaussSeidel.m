clear all;
clc
method='QUICK';
N_v=80; %Number of control volumes
rho = 1;
Gamma = 0.01;
delta = 1/N_v;
x = zeros(1,N_v+2);
for i=2:N_v+2
    if(i==2 || i==N_v+2)
        x(i)=x(i-1)+delta/2;
    else
        x(i)=x(i-1)+delta;
    end
end
    
y = x;
% boundary condition
N = N_v+2; % Number of grid points
phi_curr = zeros(N); phi_curr(:,1)=transpose(1-y);
phi_prev = zeros(N); phi_prev(:,1)=transpose(1-y);
phi_curr(N,1)=0;phi_prev(N,1)=0;
AP=zeros(N); AE=zeros(N); AW=zeros(N); AN=zeros(N);AS=zeros(N);
QP = zeros(N);
% For every node point P store A_E, A_S, A_P, A_W and A_N in matrix. 
% Also store right hand side term Qp (though it is 0 in this problem).

for i=2:N-1
    for j=2:N-1
          if(j==2) %West boundary
             AW(i,j)=min(-rho*x(j-1)*delta,0) -Gamma*2;
             AP(i,j)=AP(i,j)+Gamma*2;
          else
              AW(i,j)=min(-rho*x(j-1)*delta,0)-Gamma;
              AP(i,j)=AP(i,j)-AW(i,j);
          end
          if(i==N-1) %North boundary
              AN(i,j)=(min(-rho*y(i+1)*delta,0))-Gamma*2;
              AP(i,j)=AP(i,j)+Gamma*2;
          else
              AN(i,j)=(min(-rho*y(i+1)*delta,0))-Gamma;
              AP(i,j)=AP(i,j)-AN(i,j);
          end
          if(j==N-1) %East boundary
              AE(i,j)=0;
              AP(i,j)=AP(i,j)+max(rho*x(j+1)*delta,0);
          else
              AE(i,j)=min(rho*x(j+1)*delta,0)-Gamma;
              AP(i,j)=AP(i,j)-AE(i,j);
          end
          if(i==2) %South boundary
              AS(i,j)=0;
              AP(i,j)=AP(i,j)+max(rho*y(i-1)*delta,0);
          else
              AS(i,j)=min(rho*y(i-1)*delta,0)-Gamma;  
              AP(i,j)=AP(i,j)-AS(i,j);
          end       
    end
end


%% Solving the system of equations
iter=0;
err = 10000;
erros=[];
while iter<10000 && err > 1e-5
    for i=N-1:-1:2
        for j=2:N-1
            if(i~=2 && j~=N-1)
                phi_curr(i,j) = QP(i,j) - AW(i,j)*phi_curr(i,j-1)...
                    -AS(i,j)*phi_curr(i-1,j)-AN(i,j)*phi_curr(i+1,j) - AE(i,j)*phi_curr(i,j+1);
                phi_curr(i,j)=phi_curr(i,j)/AP(i,j);
            elseif(i~=2 && j==N-1)
                %At the east bondary
                phi_curr(i,j) = QP(i,j) - AW(i,j)*phi_curr(i,j-1)...
                    -AS(i,j)*phi_curr(i-1,j)-AN(i,j)*phi_curr(i+1,j);
                phi_curr(i,j)=phi_curr(i,j)/(AP(i,j)+AS(i,j));
                phi_curr(i,j+1)=phi_curr(i,j);
            elseif(i==2 && j~=N-1)
                %At the south boundary
                phi_curr(i,j) = QP(i,j) - AW(i,j)*phi_curr(i,j-1)...
                   -AN(i,j)*phi_curr(i+1,j) - AE(i,j)*phi_curr(i,j+1);
                phi_curr(i,j)=phi_curr(i,j)/(AP(i,j)+AS(i,j));
                phi_curr(i-1,j)=phi_curr(i,j);
            else
                %At both south and east boundary
                 phi_curr(i,j) = QP(i,j) - AW(i,j)*phi_curr(i,j-1)...
                   -AN(i,j)*phi_curr(i+1,j);
                 phi_curr(i,j)=phi_curr(i,j)/(AP(i,j)+AS(i,j)+AE(i,j));
                 phi_curr(i-1,j)=phi_curr(i,j);
                 phi_curr(i,j+1)=phi_curr(i,j);
            end
        end
    end
    err = norm(phi_prev-phi_curr);
    iter=iter+1;
    error(iter)=err;
    phi_prev = phi_curr;
end

s=contour(x,y,phi_curr,'linewidth',0.75);
title(sprintf('Contour plot of phi by %s for %iX%i',method,N_v,N_v));
% saveas(gcf, sprintf('%s_Gamma_%.3f_%iX%i.png', method,Gamma,N_v,N_v));