function [Diff_x,Diff_y,Conv_x,Conv_y]=Conv_Diff(rho,nu,delta,u_prev,v_prev,N)
% Computing for all inside nodes
Diff_x=zeros(N);Diff_y=zeros(N);Conv_x=zeros(N);Conv_y=zeros(N);

ue=(u_prev(2:N-1,2:N-1)+u_prev(2:N-1,3:N))/2; ve=(v_prev(2:N-1,2:N-1)+v_prev(2:N-1,3:N))/2;
uw=(u_prev(2:N-1,1:N-2)+u_prev(2:N-1,2:N-1))/2;vw=(v_prev(2:N-1,1:N-2)+v_prev(2:N-1,2:N-1))/2;
un=(u_prev(2:N-1,2:N-1)+u_prev(3:N,2:N-1))/2; vn=(v_prev(2:N-1,2:N-1)+v_prev(3:N,2:N-1))/2;
us=(u_prev(2:N-1,2:N-1)+u_prev(1:N-2,2:N-1))/2; vs=(v_prev(2:N-1,2:N-1)+v_prev(1:N-2,2:N-1))/2;

Conv_x(2:N-1,2:N-1)= (ue.*ue -uw.*uw + un.*vn - us.*vs)*(rho*delta);
Conv_y(2:N-1,2:N-1) = (ue.*ve-uw.*vw + vn.*vn - vs.*vs)*(rho*delta);

Diff_x(2:N-1,2:N-1)=(u_prev(2:N-1,3:N)- u_prev(2:N-1,2:N-1)...     % east face uE-uP 
                     -(u_prev(2:N-1,2:N-1)-u_prev(2:N-1,1:N-2))...   % west face uP-uW
                     +u_prev(3:N,2:N-1)-u_prev(2:N-1,2:N-1)...     % north face uN-uP
                     -(u_prev(2:N-1,2:N-1)-u_prev(1:N-2,2:N-1)) )...  % south face uP-uS
                     *(delta/delta)*nu; 
                 
Diff_y(2:N-1,2:N-1)=(v_prev(2:N-1,3:N)- v_prev(2:N-1,2:N-1)...     % east face vE-vP 
                     -(v_prev(2:N-1,2:N-1)-v_prev(2:N-1,1:N-2))...   % west face vP-vW
                     +v_prev(3:N,2:N-1)-v_prev(2:N-1,2:N-1)...     % north face vN-vP
                     -(v_prev(2:N-1,2:N-1)-v_prev(1:N-2,2:N-1)))...  % south face vP-vS
                     *(delta/delta)*(nu*rho);
end