clear all;
clc;
load('mapQUICK.mat')
load('mapObj.mat')
load('mapCDS.mat')

fluxCDS=[];
fluxUDS=[];
fluxQUICK=[];
CVs=[];
i=0;
% Loop for CDS
for k=keys(mapCDS)
    i=i+1;   
    thekey=k{1};
    CVs(i)=thekey*thekey; 
    
    x=mapCDS(thekey);
    fluxCDS(i)=2*sum(x(:,1)-x(:,2));
    x=mapObj(thekey);
    fluxUDS(i)=2*sum(x(:,1)-x(:,2));
    x=mapQUICK(thekey);
    fluxQUICK(i)=2*sum(x(:,1)-x(:,2));
end

for i=1:length(fluxCDS)
   if(abs(fluxCDS(i))>1e70)
       continue;
   else
       ite=i;
       break;
   end
end

loglog(CVs(ite:length(CVs)),fluxCDS(ite:length(fluxCDS)), 'LineWidth',1)
% hold on
ylabel('Flux');
xlabel('Number of CVs');
grid on
saveas(gcf, sprintf('Plot_of_flux_vs_CVs_CDS.png'))
loglog(CVs,fluxUDS, 'LineWidth',1)
% hold on
ylabel('Flux');
xlabel('Number of CVs');
grid on
saveas(gcf, sprintf('Plot_of_flux_vs_CVs_UDS.png'))

loglog(CVs,fluxQUICK,'LineWidth',1)
ylabel('Flux');
xlabel('Number of CVs');
grid on
saveas(gcf, sprintf('Plot_of_flux_vs_CVs_QUICK.png'))
% title(sprintf('Plot of flux of phi vs. number CVs'));
% legend('CDS', 'UDS', 'QUICK');

% saveas(gcf, sprintf('Plot_of_flux_vs_CVs.png'))