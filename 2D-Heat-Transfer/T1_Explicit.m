alpha = 0.835;
tl = 0;               %Lower bound for time interval
tu = 100;            %Upper bound for time interval
xu = 40;              %Upper bound on x; Length of plate
xl = 0;               %Lower bound of x

dt = 10;
dx = 10;              %dy = dx
lambda = (alpha*dt)/(dx^2);
t = zeros((tu-tl)/dt + 1);
%1st-dimesion is x (ROW); 2-nd is y (Column)
T_explicit = zeros((xu-xl)/dx + 1, (xu-xl)/dx + 1, (tu-tl)/dt + 1);
T_explicit(:, 1,:) = 75;                  %The left-edge
T_explicit(:, (xu-xl)/dx + 1, :) = 50; %The right-edge
T_explicit(1,:, :) = 0;                  %The bottom edge
T_explicit((xu-xl)/dx + 1, :, :) = 100;   %The top edge

%Loop to do explicite method; for every time
%Looping over time
for n=2:(tu-tl)/dt+1
    %Loop over all interior nodes
    %i is the index of x; the column no.
    for j=2:(xu-xl)/dx
        % j is y; the row number
        for i=2:(xu-xl)/dx
            T_explicit(j,i,n) = T_explicit(j,i,n-1) + lambda*(T_explicit(j,i+1,n-1) + T_explicit(j,i-1,n-1) + T_explicit(j-1,i,n-1) + T_explicit(j+1,i,n-1) - 4*T_explicit(j,i,n-1));  
        end
    end
end
%%%%%%%%%%%% TO DO PLOTTING %%%%%%%%%%%%%
% T_last = T_explicit(:,:,11);
% [x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
% [C,h] = contour(x,y,T_last);

%%%%%%%%%%%%%%%%%%%%%%%%
% T_last = T_explicit(:,:,11);
% [x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
% [C,h] = contour3(x,y,T_last,33);
% clabel(C,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% T_last = T_explicit(:,:,11);
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

T_last = T_explicit(:,:,11);
[x,y] = meshgrid(xl:dx:xu, xl:dx:xu); 
mesh(x,y,T_last);
%plottools('on')