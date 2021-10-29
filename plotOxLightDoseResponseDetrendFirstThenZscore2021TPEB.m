


%% BASELINE: takes detrended trace and performs zScore


for mouseLoop=1:size(rawDETRENDED,4)

    baseMat=[];
    
    tempMainMat=[];
    
    tempMainMat=squeeze(rawDETRENDED(:,:,:,mouseLoop));
    
    if baseMatVal==0;
    
    baseMat=squeeze(mean(tempMainMat(baseStart:baseStop,:,:)));
    
    elseif baseMatVal==1;
        
        baseMat=squeeze(median(tempMainMat(baseStart:baseStop,:,:)));
    
    elseif baseMatVal==2;
        
        baseMat=squeeze(max(tempMainMat(baseStart:baseStop,:,:)));
        
    end
    
    stdMat=squeeze(std(tempMainMat(baseStart:baseStop,:,:)));
    
    
    for sampleLoop=1:size(tempMainMat,1)
    
    tempMainMatNorm(sampleLoop,:,:)=((squeeze(tempMainMat(sampleLoop,:,:))-baseMat)./stdMat);
      
    end
    
    
    dFoverFBaselineOffsetMatDetrendFirst(:,:,:,mouseLoop)=tempMainMatNorm;

end


dFoverFMatBaselineOffsetDETRENDED=mainTrialDetrendMat;

dFoverFMatBaselineOffset_detrendFIRSTName=strcat('dFoverFBaselineOffsetMatDetrendFirst', titleStub);
save(dFoverFMatBaselineOffset_detrendFIRSTName, 'dFoverFBaselineOffsetMatDetrendFirst' );


%%

meanDoseFoverFMatBaselineOffsetDetrendFirst=squeeze(mean(dFoverFBaselineOffsetMatDetrendFirst,2)); %get trial averages

meanDoseFoverFBaselineOffsetMouseMatDetrendFirst=squeeze(mean(meanDoseFoverFMatBaselineOffsetDetrendFirst,3)); % get mouse averages


%%

if getMeanOrMaxVals==0;
    
    valMat=squeeze(mean(dFoverFBaselineOffsetMatDetrendFirst(stimValStart:stimValStop,:,:,:)));
    
elseif getMeanOrMaxVals==1;
    
    valMat=squeeze(max(dFoverFBaselineOffsetMatDetrendFirst(stimValStart:stimValStop,:,:,:)));
    
end

%%

reshapeMat=[];


if size(valMat,3) == 2;

    reshapeMat=[squeeze(valMat(:,:,1));squeeze(valMat(:,:,2))]; %squeeze(valMat(:,:,3))];

elseif size(valMat,3) == 3;
    
    reshapeMat=[squeeze(valMat(:,:,1));squeeze(valMat(:,:,2)); squeeze(valMat(:,:,3))];
    
end

%%

save('reshapeValMat_detrendFIRST', 'reshapeMat' );

%% plot traces and values

figure('Renderer', 'painters', 'Position', [1300 700 800 350]);

subplot(1,2,1)

for plotLoop=1:size(meanDoseFoverFBaselineOffsetMouseMatDetrendFirst,2);


    
    probs=plotLoop;

    plot(meanDoseFoverFBaselineOffsetMouseMatDetrendFirst(:,probs), plotColList{probs})
    hold on;
    
end

yL = get(gca,'YLim');

line([stimStart stimStart],[yL(1) yL(2)],'Color','r');


if forceAxis==1;
    
   axis([0 size(meanDoseFoverFBaselineOffsetMouseMatDetrendFirst,1) yLForceTraces(1) yLForceTraces(2)])
    
end

xlabel('zScore DETREND first')


subplot(1,2,2)

for plotLoop=1:size(meanDoseFoverFBaselineOffsetMouseMatDetrendFirst,2);
    
    
    probs=plotLoop;

    
    plot(plotLoop,reshapeMat(:,probs),plotColListVal{probs},'MarkerSize',1)
    hold on;
    
    
    
    if meanOrMedianOfMice==0;
    
    plot(plotLoop,squeeze(mean(reshapeMat(:,probs))), plotColListVal{probs},'MarkerSize',5)
    
    elseif meanOrMedianOfMice==1;
   
    plot(plotLoop,squeeze(median(reshapeMat(:,probs))), plotColListVal{probs},'MarkerSize',5)
        
    end 
    
    
    
end

yL = get(gca,'YLim');

axis([0 6 yL(1) yL(2)])

xlabel('values during stim')


if forceAxis==1;
    
   axis([0 6 yLForceVals(1) yLForceVals(2)])
    
end


if meanOrMedianOfMice==0;
    supPlotTitle=strcat('detrendedZscore', 'MEAN', titleStub);
elseif meanOrMedianOfMice==1;
    supPlotTitle=strcat('detrendedZscore', 'MEDIAN', titleStub);
end


suptitle(supPlotTitle);

saveas(gcf,supPlotTitle,'png');
saveas(gcf,supPlotTitle,'svg');

