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
