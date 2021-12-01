function P_curr = ADI_helper_naive(P,u_pred,v_pred,rho,delta_t,delta,N)
% P_curr=P;
P_curr=zeros(N+2);

% Do a column sweep
% column sweep, index i runs from row 2 to row M-1
for i=2:N-1
    
    Ap = zeros(N-2,1)-4;
    Ae = zeros(N-2-1,1)+1;
    Aw = zeros(N-2-1,1)-1;
    RHS=(2*delta*rho/(delta_t))*(u_pred(2:N-1,i+1)-u_pred(2:N-1,i-1) ...
        +v_pred(3:N,i)-v_pred(1:N-2,i));
            
    Ap(1)=Ap(1)-1;
    Ap(N-2)=Ap(N-2)+1;
    
    
    if(i==N-1)
        Ap=Ap+1;
        B=RHS+P(2:N-1,i-2);
    elseif(i==2)
        Ap=Ap-1;
        B=RHS-P(2:N-1,i+2);
    else
        B=RHS-P(2:N-1,i+2)+P(2:N-1,i-2);
    end

    P_curr(3:N,i)=transpose(ThomasAlgo(Ap,Aw,Ae,B));
    P_curr(2,i)=P_curr(3,i);
    P_curr(1,i)=P_curr(3,i);
    P_curr(N+1,i)=P_curr(N,i);
    P_curr(N+2,i)=P_curr(N,i);
    P(:,i)=P_curr(2:N+1,i);
end
P_curr(3:N,1)=P_curr(3:N,2);
P_curr(3:N,N)=P_curr(3:N,N-1);

P=P_curr;

%Do a row sweep
% row sweep, index i runs from row 2 to row N-1
for j=N-1:-1:2
    Ap = zeros(N-2,1)-4;
    Ae = zeros(N-2-1,1)+1;
    Aw = zeros(N-2-1,1)-1;
    RHS=(2*delta*rho/(delta_t))*(u_pred(j,3:N)-u_pred(j,1:N-2) ...
         +v_pred(j+1,2:N-1)-v_pred(j-1,2:N-1));
            
    Ap(1)=Ap(1)+1;
    Ap(N-2)=Ap(N-2)-1;

    if(j==N-1)
        Ap=Ap+1;
        B=RHS+P(j-1,2:N-1);
    elseif(j==2)
        Ap=Ap+1;
        B=RHS-P(j+1,2:N-1);
    else
        B=RHS-P(j+1,2:N-1)+P(j-1,2:N-1);
    end
    
    P_curr(j,3:N)=transpose(ThomasAlgo(Ap,Aw,Ae,B));
    P_curr(j,2)=P_curr(j,3);
    P_curr(j,1)=P_curr(j,3);
    P_curr(j,N+1)=P_curr(j,N);
    P_curr(j,N+2)=P_curr(j,N);
    P(j,:)=P_curr(j,:);
end
P_curr(1,3:N)=P_curr(2,3:N);
P_curr(N,3:N)=P_curr(N-1,3:N);

P_curr=P_curr(2:N+1,2:N+1);

end