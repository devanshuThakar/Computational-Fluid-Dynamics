alpha = 0.835;
L = 10; %m

tl = 0;
tu = 100;
xu = 40;
xl = 0;

dt = 10;
dx = 10; %dy = dx
lambda = (alpha*dt)/(dx^2);
t = zeros((tu-tl)/dt + 1);
%1st-dimesion is x (ROW); 2-nd is y (Column); 3rd is time t
T = zeros((xu-xl)/dx + 1, (xu-xl)/dx + 1, (tu-tl)/dt + 1);
T(:, 1,:) = 75;                     %The left-edge
T(:, (xu-xl)/dx + 1, :) = 50;       %The right-edge
T(1,:, :) = 0;                    %The top edge
T((xu-xl)/dx + 1, :, :) = 100;        %The Bottom edge

%Loop to do the Column sweep version of ADI method
%Looping over time
%To store the value of Temperature at n+1/2, a temperory array is created
for n=2:(tu-tl)/dt+1
    % i,j in report represents in code, ith-column,j-row
    T_half = T(:, :, n-1);
    %Iterate over all i(columns) as in equation-1
    for i=2:(xu-xl)/dx
        R = zeros((xu-xl)/dx - 1, 1);
        f = zeros((xu-xl)/dx - 1, 1);
        g = zeros((xu-xl)/dx - 2, 1);
        e = zeros((xu-xl)/dx - 2, 1);
        f = f + 2*(1+lambda);
        g = g - lambda;
        e = e - lambda;
        for j=2:(xu-xl)/dx
            R(j-1) = lambda*T_half(j,i-1) + 2*(1-lambda)*T_half(j,i) + lambda*T_half(j,i+1);
        end
        R(1) = R(1) + lambda*T(1,i,n-1);
        R((xu-xl)/dx - 1) = R((xu-xl)/dx - 1) + lambda*T((xu-xl)/dx+1,i,n-1);
        T_half(2:(xu-xl)/dx,i) = ThomasAlgorithm(e,f,g,R);
    end
    T(:,:,n) = T_half;
    %Iterate over all j(rows) as in equation-2
    for j=(xu-xl)/dx:-1:2
        for i=2:(xu-xl)/dx
            R(i-1) = lambda*(T(j-1,i,n)) + (2*(1-lambda)*T(j,i,n)) + lambda*(T(j+1,i,n));
        end
        R(1) = R(1) + (lambda*T_half(j,1));
        R((xu-xl)/dx - 1) = R((xu-xl)/dx - 1) + (lambda*T_half(j,(xu-xl)/dx + 1));
        T(j,2:(xu-xl)/dx,n) = transpose(ThomasAlgorithm(e,f,g,R));
    end
end

%%%%%%%%%%%% TO DO PLOTTING %%%%%%%%%%%%%
% T_last = T(:,:,(tu-tl)/dt + 1);
% [x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
% [C,h] = contour(x,y,T_last,20);

% v = VideoWriter('Implicit_Counter_Fine.avi');
% v.FrameRate = 3;
% open(v);
% [x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
% for k=1:(tu-tl)/dt + 1
%    [C,h] = contourf(x,y,T(:,:,k));
%    frame = getframe(gcf);
%    writeVideo(v,frame);
%    pause(0.15)
% end
% close(v);

%%%%%%%%%%%%%%%%%%%%%%%%
% T_last = T(:,:,(tu-tl)/dt+1);
% [x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
% [C,h] = contour3(x,y,T_last);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% T_last = T(:,:,(tu-tl)/dt + 1);
% [x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
% [C,h] = contourf(x,y,T_last);
% clabel(C,h)
%%%%%%%% Video %%%%%%%%%%%%%%
% v = VideoWriter('exp3Video.avi');
% v.FrameRate = 3;
% open(v);
% [x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
% for k=1:(tu-tl)/dt + 1
%    [C,h] = contourf(x,y,T_explicit(:,:,k));
%    frame = getframe(gcf);
%    writeVideo(v,frame);
%    pause(0.15)
% end
% close(v);

T_last = T(:,:,(tu-tl)/dt  + 1);
[x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
mesh(x,y,T_last);

%f is the digonal, array. e is the lower band digonal and g is upper one
function Thomas = ThomasAlgorithm(e,f,g,R)
    n = length(f);
    X = zeros(1, n);
    %Decomposition
    for i=2:n
        %index of e,will be i-1
        e(i-1) = e(i-1)/f(i-1);
        %index of g,will be i
        f(i) = f(i) - e(i-1)*g(i-1);
    end
    %Forward substitution
    for i=2:n
        R(i) = R(i) - R(i-1)*e(i-1);
    end
    % Backward Substitution to get X
    X(n) = R(n)/f(n);
    for i=n-1:-1:1
        X(i) = (R(i) - g(i)*X(i+1))/f(i);
    end
    Thomas = X;
end

