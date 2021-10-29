

load('mainOxCtrlRunning.mat');
load('mainOxCtrlRunningCtrl.mat');

%%

mainRunningMat=mainOxCtrlRunning;
mainRunningMatCtrl=mainOxCtrlRunningCtrl;

%%

runStart=125;
runStop=175;

meanOrMax=0; %0 = max, 1 = mean 2= median

forceAxis=1;

yRange=[-5 60];
xRange=[0 1];

%%

if meanOrMax==0;


maxRun=squeeze(max(mainRunningMat(runStart:runStop,1,:,:)));
maxOrx=squeeze(max(mainRunningMat(runStart:runStop,2,:,:)));


elseif meanOrMax==1;

maxRun=squeeze(mean(mainRunningMat(runStart:runStop,1,:,:)));
maxOrx=squeeze(mean(mainRunningMat(runStart:runStop,2,:,:)));


elseif meanOrMax==2;

maxRun=squeeze(median(mainRunningMat(runStart:runStop,1,:,:)));
maxOrx=squeeze(median(mainRunningMat(runStart:runStop,2,:,:)));


elseif meanOrMax==3;

maxRun=squeeze(mainRunningMat(runStart:runStop,1,:,:));
maxOrx=squeeze(mainRunningMat(runStart:runStop,2,:,:));


end


%%

if meanOrMax==0;


maxRunCtrl=squeeze(max(mainRunningMatCtrl(runStart:runStop,1,:,:)));
maxOrxCtrl=squeeze(max(mainRunningMatCtrl(runStart:runStop,2,:,:)));


elseif meanOrMax==1;

maxRunCtrl=squeeze(mean(mainRunningMatCtrl(runStart:runStop,1,:,:)));
maxOrxCtrl=squeeze(mean(mainRunningMatCtrl(runStart:runStop,2,:,:)));


elseif meanOrMax==2;

maxRunCtrl=squeeze(median(mainRunningMatCtrl(runStart:runStop,1,:,:)));
maxOrxCtrl=squeeze(median(mainRunningMatCtrl(runStart:runStop,2,:,:)));


elseif meanOrMax==3;

maxRunCtrl=squeeze(mainRunningMatCtrl(runStart:runStop,1,:,:));
maxOrxCtrl=squeeze(mainRunningMatCtrl(runStart:runStop,2,:,:));


end


%%

maxOrx=maxOrx(:);
maxRun=maxRun(:);

maxOrxCtrl=maxOrxCtrl(:);
maxRunCtrl=maxRunCtrl(:);


%%

figure;

[rVal,pVal]=corrcoef(maxRun,maxOrx);
rVal=rVal.^2;

pFit = polyfit(maxRun,maxOrx,1);

f = polyval(pFit,maxRun);

plot(maxRun,maxOrx,'og',maxRun,f,'-k','markerSize',3);
hold on

rOxRun=rVal(2);
pOxRun=pVal(2);

txt=strcat('r2 OxRun = ', num2str(rOxRun), ' p OxRun = ', num2str(pOxRun));
 text(0, 5, txt);


[rVal,pVal]=corrcoef(maxRunCtrl,maxOrxCtrl);
rVal=rVal.^2;

pFit = polyfit(maxRunCtrl,maxOrxCtrl,1);

f = polyval(pFit,maxRun);

plot(maxRunCtrl,maxOrxCtrl,'ok',maxRun,f,'-k','markerSize',3);
hold on

rOxRunCtrl=rVal(2);
pOxRunCtrl=pVal(2);

txt=strcat('r2 OxRunCtrl = ', num2str(rOxRunCtrl), ' p OxRunCtrl = ', num2str(pOxRunCtrl));
 text(0, 7, txt); 
 

 if forceAxis==1
     
     axis([xRange(1) xRange(2) yRange(1) yRange(2)])
     
 else        
 end
 
 ylabel('Max Z-Score(SD)')
 
 
 
 if meanOrMax==0;
 
     xlabel(strcat('running MAX',num2str(runStart),'to', num2str(runStop)))
     
 elseif meanOrMax==1;
     
     xlabel(strcat('running MEAN',num2str(runStart),'to', num2str(runStop)))
     
 elseif meanOrMax==2;
 
     xlabel(strcat('running MEDIAN',num2str(runStart),'to', num2str(runStop)))
     
 end

 
saveas(gcf,'corr oxLight vs Run','svg');
saveas(gcf,'corr oxLight vs Run','png');


