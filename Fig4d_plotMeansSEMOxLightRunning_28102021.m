
for runLoop=1:2


    if runLoop==1;
    
        mainOxCtrlRunningMatLoad=load('mainOxCtrlRunning.mat');
        mainOxCtrlRunningMat=mainOxCtrlRunningMatLoad.mainOxCtrlRunning;
        
        titleStub='oxLight';
        graphShiftVal=0;
            
    elseif runLoop==2;
        
        mainOxCtrlRunningMatLoad=load('mainOxCtrlRunningCtrl.mat');
        mainOxCtrlRunningMat=mainOxCtrlRunningMatLoad.mainOxCtrlRunningCtrl;
        titleStub='oxLightControl';
        graphShiftVal=300;
        
    end

%%

runOnset=125;

meanRunningOfRepeats=squeeze(mean(mainOxCtrlRunningMat,3));

meanRunningRepeatsMice=squeeze(mean(meanRunningOfRepeats,3));

errRunningRepeatsMice=std(meanRunningOfRepeats,[],3)/sqrt(size(meanRunningOfRepeats,3)); 

%% values for Graph

forceAxis=1;
yRunMin=0;
yRunMax=0.16;

yOxMin=-1;
yOxMax=8;

plotRange=75:220; 

%% plots mean +/- SEM traces




figure('Renderer', 'painters', 'Position', [150+graphShiftVal 150 300 700]);

subplot(2,1,1)

channel=1;
shadedErrorBar(plotRange,meanRunningRepeatsMice(plotRange,channel),errRunningRepeatsMice(plotRange,channel),'lineprops','-m','transparent',1)
hold on;
plot(plotRange,meanRunningRepeatsMice(plotRange,channel),'-k')

ylabel('running')

if forceAxis==0

yL = get(gca,'YLim');

elseif forceAxis ==1;
    
yL = [yRunMin, yRunMax];

end

line([runOnset runOnset],[yL(1) yL(2)],'Color','b');
axis([plotRange(1) plotRange(end) yL(1) yL(2)])


subplot(2,1,2)

channel=2;
shadedErrorBar(plotRange,meanRunningRepeatsMice(plotRange,channel),errRunningRepeatsMice(plotRange,channel),'lineprops','-g','transparent',1)
hold on;
plot(plotRange,meanRunningRepeatsMice(plotRange,channel),'-k')

if forceAxis==0

yL = get(gca,'YLim');

elseif forceAxis ==1;
    
yL = [yOxMin, yOxMax];

end

line([runOnset runOnset],[yL(1) yL(2)],'Color','b');
axis([plotRange(1) plotRange(end) yL(1) yL(2)])

ylabel('oxLight')

SaveSupTitle=strcat('runOxLightmeanSEM',titleStub);

suptitle(SaveSupTitle);

saveas(gcf,SaveSupTitle,'png');
saveas(gcf,SaveSupTitle,'svg');



end
