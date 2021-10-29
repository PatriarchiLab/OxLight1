

load('mainOxCtrlRunning.mat');
load('mainOxCtrlRunningCtrl.mat');
        
for channelLoop=1:2

channel = channelLoop;

if channel == 1 
    
titleStub = 'running speed degSec';   
graphShiftVal=0; 


elseif channel == 2
    
titleStub = 'oxLight MaxZscoreSD';
graphShiftVal=300;

end

startVal=125;
stopVal=175;

oxLightCtrlMat=squeeze(max(mainOxCtrlRunningCtrl(startVal:stopVal,channel,:,:)));

oxLightMat=squeeze(max(mainOxCtrlRunning(startVal:stopVal,channel,:,:)));


%%

oxLightVect=oxLightMat(:);
oxLightCtrlVect=oxLightCtrlMat(:);

%%

meanOxVal=mean(oxLightVect);
errOxVal=squeeze(nanstd(oxLightVect,[],1)/sqrt(size(oxLightVect,1)));


meanOxValCtrl=mean(oxLightCtrlVect);
errOxValCtrl=squeeze(nanstd(oxLightCtrlVect,[],1)/sqrt(size(oxLightCtrlVect,1)));

%% find rank sum p values 

[p r]=ranksum(oxLightVect,oxLightCtrlVect);

%% plots single values, and mean +/- SEM

figure('Renderer', 'painters', 'Position', [350+graphShiftVal 350 300 600]);

plot(1,oxLightVect,'.g')
hold on
errorbar(0.75,meanOxVal,errOxVal,'og')

plot(2,oxLightCtrlVect,'.k')
errorbar(2.5,meanOxValCtrl,errOxValCtrl,'ok')


yL = get(gca,'YLim');

%axis([0 3 yL(1) yL(2)])

    longbar = {'oxLight','oxLightCtrl'};
    set(gca,'xtick',1:2,'xticklabel',longbar)

ylabel(titleStub)

suptitleText=strcat('means SEM, and individual trials',titleStub);

suptitle('means SEM, and individual trials')





end
