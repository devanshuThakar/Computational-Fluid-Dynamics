% Function to solve linear algebraic equation system using Gauss elimination:
% Here, A: coefficient matrix
%       b: RHS vector
%
%Author: Kausadikar Varad Prashant


function Y = Gauss_elim_fun(A,b)
    
    n=length(b);
% Augmented matrix:
    Ab=[A, b];
%% Forward elimination:

for k=1:(n-1)
    
    for i=(k+1):n  
        
        Ab(i,k:n+1)=Ab(i,k:n+1)-(Ab(i,k)/Ab(k,k))*Ab(k,k:n+1);

    end
end

%% Backward substitution:

    %Solution vector:
        Y=NaN(n,1);
 

for i=n:-1:1
      
   Y(i,1)=( Ab(i,end)-Ab(i,i+1:n)*Y(i+1:n,1) ) / Ab(i,i);
        
end
    
end

