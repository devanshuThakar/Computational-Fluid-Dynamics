clear all;
clc;
method='Rhie_chow';
Nv=80;             %Number of control volumes
delta = 1/Nv;
nu=0.01;           %kinematic viscosity
rho=1;
t_start=0;
t_end=15;
t=t_start;
N_t=1000;          %Number of time steps
% delta_t=(t_end-t_start)/N_t;
dt=0.001;

steady_state=false;

x=zeros(1,Nv+2);
for i=2:Nv+2
    x(i)=x(i-1)+delta;
end
y = x;
N = Nv+2; % Number of grid points including ghost cells

u_curr=zeros(N);v_curr=zeros(N);P_curr=zeros(N);
% Setting the top boundary condition
u_curr(N,:)=2;u_curr(:,1)=0;u_curr(:,N)=0;
u_prev=u_curr;v_prev=v_curr;P_prev=P_curr;

% Intermediate or predicted velocity
v_pred=v_prev;u_pred=u_prev;
% A varable to hold velocity temporary
u_temp=zeros(N);
v_temp=zeros(N);

% residual=[];
time=[];
u_(:,:,1)=u_curr;
v_(:,:,1)=v_curr;
time(1)=t;
cnt=0;

while t<t_end && ~steady_state
 % Assigning predicted equal to previous only for programmatic simplicity
 u_pred=u_prev;
 v_pred=v_prev;
 P_pred=P_prev;
 
 converge_Rhi_chow=false;
 iter=0;
 err_prev=-10;
 
 while(~converge_Rhi_chow && iter<5000)
     
 [Diff_x,Diff_y,Conv_x,Conv_y]=Conv_Diff(rho, nu,delta,u_pred,v_pred,N);
  
 u_pred(2:N-1,2:N-1) =  u_prev(2:N-1,2:N-1)+ (Diff_x(2:N-1,2:N-1)-Conv_x(2:N-1,2:N-1))*(dt/(rho*delta*delta));
 v_pred(2:N-1,2:N-1) = v_prev(2:N-1,2:N-1) + (Diff_y(2:N-1,2:N-1)-Conv_y(2:N-1,2:N-1))*(dt/(rho*delta*delta));
 
  % Update boundary conditions for pred
 u_pred(N,2:N-1)=2*1-u_pred(N-1,2:N-1); v_pred(N,2:N-1)=(2*0)-v_pred(2,2:N-1);
 u_pred(1,2:N-1)=(2*0)-u_pred(2,2:N-1); v_pred(1,2:N-1)=(2*0)-v_pred(2,2:N-1);
 u_pred(2:N-1,N)=(2*0)-u_pred(2:N-1,N-1);v_pred(2:N-1,N)=(2*0)-v_pred(2:N-1,N-1);
 u_pred(2:N-1,1)=(2*0)-u_pred(2:N-1,2); v_pred(2:N-1,1)=(2*0)-v_pred(2:N-1,2);

 P_pred=ADI_helper(P_pred,u_pred,v_pred,rho,dt,delta,N);
 
 %Velocity update
 u_curr(2:N-1,2:N-1)=u_pred(2:N-1,2:N-1)- (dt/(rho*delta))*(P_pred(2:N-1,3:N)-P_pred(2:N-1,2:N-1));
 v_curr(2:N-1,2:N-1)=v_pred(2:N-1,2:N-1)- (dt/(rho*delta))*(P_pred(3:N,2:N-1)-P_pred(2:N-1,2:N-1));

 % Update boundary conditions for current velocity
 u_curr(N,2:N-1)=(2*1)-u_curr(N-1,2:N-1); v_curr(N,2:N-1)=(2*0)-v_curr(2,2:N-1);
 u_curr(1,2:N-1)=(2*0)-u_curr(2,2:N-1); v_curr(1,2:N-1)=(2*0)-v_curr(2,2:N-1);
 u_curr(2:N-1,N)=(2*0)-u_curr(2:N-1,N-1);v_curr(2:N-1,N)=(2*0)-v_curr(2:N-1,N-1);
 u_curr(2:N-1,1)=(2*0)-u_curr(2:N-1,2); v_curr(2:N-1,1)=(2*0)-v_curr(2:N-1,2);
 
 err=norm(u_curr-u_pred)+norm(v_curr-v_pred);
 if(err<0.1 || norm((u_curr-u_temp) + (v_curr-v_temp))<1e-03 || abs(err)>abs(err_prev) || abs(abs(err)-abs(err_prev))<1e-08)
     converge_Rhi_chow=true;
 end 
 err_prev=err;
 
 u_pred=u_curr;
 v_pred=v_curr;
 P_curr=P_pred;
 
 u_temp=u_curr;
 v_temp=v_curr;
 
 iter=iter+1;
 
 end
 
 resi=norm(v_curr-v_prev)+norm(u_curr-v_prev);
%  residual(length(residual)+1)=resi;
 if(resi<1e-03)
     steadt_state=true;
 end
 v_prev=v_curr;
 u_prev=u_curr;
 P_prev=P_curr;
 
 t=t+dt;
 
 if(mod(floor((t/dt)*0.4),10)==0)
     cnt=cnt+1;
     u_(:,:,cnt)=u_curr;
     v_(:,:,cnt)=v_curr;
     time(cnt)=t;
 end
 
end

%Changing the ghost cell values to values at boundary for plotting
x(1)=0;y(1)=0;x(N)=1;y(N)=1;
u_curr(N,:)=1;u_curr(1,:)=0;u_curr(:,1)=0;u_curr(:,N)=0;
v_curr(N,:)=0;v_curr(1,:)=0;v_curr(:,1)=0;v_curr(:,N)=0; 