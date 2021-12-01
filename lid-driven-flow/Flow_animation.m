clear all;clc;
load('Field_80X80.mat');
Field=Field_80X80;
Nv=80;
N=Nv+2;

visRate = 4; % downsample rate of the data for quiver
recordGIF = true; % set to true if you wanna make GIF
recordRate = 20;
filename = 'animation_sample.gif'; % Specify the output file name

figure
% delta=1/Nv;
% x=zeros(1,Nv+2);
% for i=2:Nv+2
%     x(i)=x(i-1)+delta;
% end
x=Field('x');y=Field('y');
% y = x;
[X,Y] = meshgrid(x,y); % cell center grid
uce=zeros(Nv+2);vce=zeros(Nv+2);
[~,h_abs] = contourf(X',Y',sqrt(uce.^2+vce.^2)); % contour

hold on
% Downsample the data（d = downsampled）
xd = x(1:visRate:Nv+2);
yd = y(1:visRate:Nv+2);
[Xd,Yd] = meshgrid(xd, yd);

uced = uce(1:visRate:end,1:visRate:end);
vced = vce(1:visRate:end,1:visRate:end);
h_quiver = quiver(Xd',Yd',uced,vced,3,'Color',[1,1,1]);

hold off
xlim([0 1]); ylim([0 1]);
harrow = annotation('textarrow',[0.3 0.7],[0.96 0.96],"LineWidth",2);

haxes = gca;
haxes.XTick = [];
haxes.YTick = [];


%% Start the simulation

initialFrame = true;

for ii = 1:length(Field('time'))
    bctop = 1; % top velocity
    
    % Update the plot at every recordRate steps
    if mod(ii,recordRate) == 0
        % get velocity at the cell center (for visualization)
        tmp=Field('u');
        uce = (tmp(:,:,ii));
        %set the boundaries
        uce(:,N)=1;uce(:,1)=0;uce(N,:)=0;uce(1,:)=0;
        uce=transpose(uce);
        tmp=Field('v');
        vce = (tmp(:,:,ii)); % v at cell center
        vce(:,N)=0;vce(:,1)=0;vce(N,:)=0;vce(1,:)=0;
        vce=transpose(vce);
        % update plot (downsample)
        h_quiver.UData = uce(1:visRate:end,1:visRate:end);
        h_quiver.VData = vce(1:visRate:end,1:visRate:end);
        h_abs.ZData = sqrt(uce.^2+vce.^2);
        
        drawnow
        
        if recordGIF
            frame = getframe(gcf);
            tmp = frame2im(frame); 
            [A,map] = rgb2ind(tmp,256); 
            if initialFrame
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
                initialFrame = false;
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);% 画像をアペンド
            end
        end
        
        pause(0.1);
    end
end
