%get Time to Peak

load('mainOxCtrlRunning.mat');
mainRunningMat=mainOxCtrlRunning;

%%

runMatPeak=squeeze(mainRunningMat(125:175,1,:,:));

[maxRunMatPeakMatVal maxRunMatPeakMatIdx]=max(runMatPeak);

runPeaksVect=maxRunMatPeakMatIdx(:)/10;

meanRunPeak=squeeze(mean(runPeaksVect));

errRunPeak=std(runPeaksVect,[],1)/sqrt(size(runPeaksVect,1)); 


runMatOrxPeak=squeeze(mainRunningMat(125:175,2,:,:));

[maxOrxRunMatPeakMatVal maxOrxRunMatPeakMatIdx]=max(runMatOrxPeak);

runPeaksOrxVect=maxOrxRunMatPeakMatIdx(:)/10;


meanOrxPeak=squeeze(mean(runPeaksOrxVect));

errOrxPeak=std(runPeaksOrxVect,[],1)/sqrt(size(runPeaksOrxVect,1));


%%

figure

plot(1, runPeaksVect,'.k')
hold on
errorbar(1.25,meanRunPeak,errRunPeak,'ok')

plot(3, runPeaksOrxVect,'.g')
errorbar(2.75,meanOrxPeak,errOrxPeak,'og')

axis([0 4 0 6])

[p r]=ranksum(runPeaksVect,runPeaksOrxVect);

txt=strcat('p  ',num2str(p))

text(2,3,txt)

xlabel('Time to Peak (Sec)')

saveas(gcf,'timeToPeakOrxVsRun','svg');
saveas(gcf,'timeToPeakOrxVsRun','png');

%%

meanRunPeak
errRunPeak

meanOrxPeak
errOrxPeak


