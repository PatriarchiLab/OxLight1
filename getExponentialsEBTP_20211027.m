
clear

%% Decide To Do Frequency Dose response or 10Hz stim train duration

hzOrTrain=0; %0 = Hz dose response, 1 = Train.

%%


%endVal=150;
stimStart=300; %defines start window for finding rise times (onsets)
slopeDown=320;
delayVal=10; %20

%sgolay  filtering parameters
linearFitType=1; % 0 takes from baseline
order = 1;
framelen = 41; %5


fileLoadName=dir('dFoverFMatBaselineOffset_detrendFIRST.mat');
fileLoadName=fileLoadName.name;
load(fileLoadName);


%% variables you cannot play with

numSecs=10.7363;

mousLabelStruct={'mouse 1','mouse 2','mouse 3'};
colourStruct={'.k','.r','.m','.c','.b','.g'};
colourStructCirc={'ok','or','om','oc','ob','og'};

offsetNum=10; %for subprocess tauDecay
medianRange=offsetNum-1;

%%

if hzOrTrain==0;

channelList={'catch', '1Hz', '5Hz','10Hz','20Hz'};

offValVect=[600,600,611.13,622.09,647.11]; % for Hz 

onsetFindOffsetVal=400; %400  %finds offset of window for detecting slope of onset


plotStartFrom=3; % stim type to calculate/plot from.


elseif hzOrTrain==1;

channelList={'1s', '5s', '10s','15','20s','30s'};

offValVect=[300+numSecs*1,300+numSecs*5,300+numSecs*10,300+numSecs*15,300+numSecs*20,300+numSecs*30]; % for 10Hz trains

onsetFindOffsetVal=[400,300+numSecs*5,300+numSecs*10,300+numSecs*15,300+numSecs*20,300+numSecs*30];

plotStartFrom=1; % stim type to calculate/plot from.


end



%% Finds Rise Times


for trialNum=1:size(dFoverFBaselineOffsetMatDetrendFirst,2);  
    
    for channelLoop=plotStartFrom:size(dFoverFBaselineOffsetMatDetrendFirst,3)
        
        figure('Renderer', 'painters', 'Position', [1+(channelLoop*200) 150 200 450]);
        
        channel=channelLoop;
        
        numMice=size(dFoverFBaselineOffsetMatDetrendFirst,4);
        
        for mouseLoop=1:numMice
            
            subplot(numMice,1,mouseLoop)
            
            %[val idx]=max(dFoverFBaselineOffsetMatDetrendFirst(:,trialNum,channel,mouseLoop));
            
            idx=onsetFindOffsetVal; %offValVectON(channelLoop)
            
            exampleTrace=[];
            
            exampleTrace=dFoverFBaselineOffsetMatDetrendFirst(stimStart:idx,trialNum,channel,mouseLoop);
            
            exampleTrace = sgolayfilt(exampleTrace,order,framelen);
            
            tauDecayEquationSolver
            
            t1VectRiseTime(channelLoop,trialNum,mouseLoop)=t1;
            
            
            tauValY=max((yy)/100)*63.2;
            
            [~,tauValX] = (min(abs(yy - tauValY)));  %tauValY
            tauVectMathildePercMaxONSET(trialNum,channelLoop,mouseLoop)=tauValX;
            
            
            
            if ~isempty(t1)
            
                yL = get(gca,'YLim');

            line([tauValX tauValX],[yL(1) yL(2)],'Color','r');    
            
                
            end
                
            
            ylabel(mousLabelStruct{mouseLoop})
            
        end
        
        xlabelTitle=strcat('stim', channelList{channel},'  ONSET trial',num2str(trialNum));
        
        xlabel(xlabelTitle)
        
               saveas(gcf,'xlabelTitle','svg');
                saveas(gcf,'xlabelTitle','png');
        
    end
    
end


%

tauVectMathildePercMaxONSETReshape=[tauVectMathildePercMaxONSET(:,:,1);tauVectMathildePercMaxONSET(:,:,2);tauVectMathildePercMaxONSET(:,:,3)]/10;


%%


for nanLoop=1:size(tauVectMathildePercMaxONSETReshape,2)


tempNaNMat=tauVectMathildePercMaxONSETReshape(:,nanLoop);
    
tempNaNIdx=find(tempNaNMat<0.5);

tempNaNMat(tempNaNIdx)=nan;

tauVectMathildePercMaxONSETReshape(:,nanLoop)=tempNaNMat;
    
end

save('tauVectMathildePercMaxONSETReshapeMat','tauVectMathildePercMaxONSETReshape')


%%

meantauVectMatLabONSETMathilde=nanmean(tauVectMathildePercMaxONSETReshape);

errTauVectMatLabONSET=squeeze(nanstd(tauVectMathildePercMaxONSETReshape,[],1)/sqrt(size(tauVectMathildePercMaxONSETReshape,1)));

%% plots rise Time values (MG)


figure;


for stimLoop=plotStartFrom:size(tauVectMathildePercMaxONSETReshape,2);

stim=stimLoop;
plot(stimLoop,tauVectMathildePercMaxONSETReshape(:,stim),colourStruct{stimLoop})
hold on;
errorbar(stimLoop+0.25,meantauVectMatLabONSETMathilde(:,stim),errTauVectMatLabONSET(:,stim),colourStructCirc{stimLoop})

end

numStimsToPlot=size(plotStartFrom:size(tauVectMathildePercMaxONSETReshape,2),2);


    longbar = channelList(plotStartFrom:numStimsToPlot+plotStartFrom-1);
    set(gca,'xtick',plotStartFrom:numStimsToPlot+plotStartFrom,'xticklabel',longbar)

xlabel('Rise Time (s)')
ylabel('Stim Type')

saveas(gcf,'matlabTauRiseTimeVals','svg');
saveas(gcf,'matlabTauRiseTimeVals','png');


%% OFFSET MATHILDE

for trialNum=1:size(dFoverFBaselineOffsetMatDetrendFirst,2);  %%MAKE LOOP
    
    for channelLoop=plotStartFrom:size(dFoverFBaselineOffsetMatDetrendFirst,3)
        
        figure('Renderer', 'painters', 'Position', [1+(channelLoop*200) 150 200 450]);
        
        channel=channelLoop;
        
        numMice=size(dFoverFBaselineOffsetMatDetrendFirst,4);
        
        for mouseLoop=1:numMice
            
            subplot(numMice,1,mouseLoop)
            
            idx=offValVect(channelLoop);
            
            [val idxMin]=min(dFoverFBaselineOffsetMatDetrendFirst(idx:end,trialNum,channel,mouseLoop));
            
            idxMin=idxMin+idx;
            
            exampleTrace=[];
            
            %gets range using fixed timing.
            %exampleTrace=dFoverFBaselineOffsetMatDetrendFirst(offValVect(channelLoop)+delayVal:offValVect(channelLoop)+slopeDown+delayVal,trialNum,channel,mouseLoop);
            
            %gets range using min
            exampleTrace=dFoverFBaselineOffsetMatDetrendFirst(offValVect(channelLoop)+delayVal:idxMin-delayVal,trialNum,channel,mouseLoop);
                        
            
            
            
            exampleTrace = sgolayfilt(exampleTrace,order,framelen);
            
            
            tauDecayEquationSolver
            
            t1VectMathilde(trialNum,channelLoop,mouseLoop)=t1;
     
            tauValY=max((yy)/100)*36.8;
            
            [~,tauValX] = (min(abs(yy - tauValY)));  %tauValY
            tauVectMathildePercMaxOFFSET(trialNum,channelLoop,mouseLoop)=tauValX;
      
            
            if ~isempty(t1)
            yL = get(gca,'YLim');

            line([tauValX tauValX],[yL(1) yL(2)],'Color','r');
            
            %line([t1 t1],[yL(1) yL(2)],'Color','r');    
            
                
            end
                
            
            ylabel(mousLabelStruct{mouseLoop})
            
        end
        
        xlabelTitle=strcat('stim',channelList{channel},'  OFFSET trial',num2str(trialNum));
        
        xlabel(xlabelTitle)
        
        saveas(gcf,'xlabelTitle','svg');
        saveas(gcf,'xlabelTitle','png');

        
    end
    
end


tauVectMathildePercMaxOFFSETReshape=[tauVectMathildePercMaxOFFSET(:,:,1);tauVectMathildePercMaxOFFSET(:,:,2);tauVectMathildePercMaxOFFSET(:,:,3)]/10;


%% remove values lower than 0.1


for nanLoop=1:size(tauVectMathildePercMaxOFFSETReshape,2)


tempNaNMat=tauVectMathildePercMaxOFFSETReshape(:,nanLoop);
    
tempNaNIdx=find(tempNaNMat<0.5);

tempNaNMat(tempNaNIdx)=nan;

tauVectMathildePercMaxOFFSETReshape(:,nanLoop)=tempNaNMat;
    
end

%%

save('tauVectMathildePercMaxOFFSETReshapeMat','tauVectMathildePercMaxOFFSETReshape')


%%

meanDFoverFBaselineOffsetMatDetrendFirst=squeeze(nanmean(dFoverFBaselineOffsetMatDetrendFirst,2));

%%

meantauVectMatLabOFFSETMathilde=nanmean(tauVectMathildePercMaxOFFSETReshape);

errTauVectMatLabOFFSET=squeeze(nanstd(tauVectMathildePercMaxOFFSETReshape,[],1)/sqrt(size(tauVectMathildePercMaxOFFSETReshape,1)));


%%

figure;


for stimLoop=plotStartFrom:size(tauVectMathildePercMaxOFFSETReshape,2);

stim=stimLoop;


plot(stimLoop,tauVectMathildePercMaxOFFSETReshape(:,stim),colourStruct{stimLoop})
hold on;
errorbar(stimLoop+0.25,meantauVectMatLabOFFSETMathilde(:,stim),errTauVectMatLabOFFSET(:,stim),colourStructCirc{stimLoop})

end

numStimsToPlot=size(plotStartFrom:size(tauVectMathildePercMaxONSETReshape,2),2);


    longbar = channelList(plotStartFrom:numStimsToPlot+plotStartFrom-1);
    set(gca,'xtick',plotStartFrom:numStimsToPlot+plotStartFrom,'xticklabel',longbar)

xlabel('Decay Time (s)')
ylabel('Stim Type')

saveas(gcf,'matlabTauDecayTimeVals','svg');
saveas(gcf,'matlabTauDecayTimeVals','png');


% %%  does matlab exponents exp1 or 2 for OFFSET
% %
% 
% for trialNum=1:size(dFoverFBaselineOffsetMatDetrendFirst,2);
%     
%     for channelLoop=1:size(dFoverFBaselineOffsetMatDetrendFirst,3)
%         
%         figure('Renderer', 'painters', 'Position', [1+(channelLoop*280) 150 250 450]);
%         
%         channel=channelLoop;
%         
%         numMice=size(dFoverFBaselineOffsetMatDetrendFirst,4);
%         
%         for mouseLoop=1:numMice
%             
%             subplot(numMice,1,mouseLoop)
%             
%             %[val idx]=max(dFoverFBaselineOffsetMatDetrendFirst(:,trialNum,channel,mouseLoop));
%             
%             idx=offValVect(channelLoop); % gives from offset of stim.
%             
%             exampleTrace=[];
%             
%             exampleTrace=dFoverFBaselineOffsetMatDetrendFirst(idx+delayVal:idx+endVal,trialNum,channel,mouseLoop);
%             
%             x=1:1:size(exampleTrace,1);
%             
%             x=x';
%             
%             f=fit(x,exampleTrace,'exp2'); %exp1
%             
%             coeffVals=coeffvalues(f)
%             tauValXCoeff=coeffVals(2); % gets from matlab calculated values.
%             
%             yVals = f(x);
%             tauValY=max((yVals)/100)*36.8;
%             
%             [~,tauValX] = (min(abs(yVals - tauValY)));  %tauValY
%             tauVectMatLabOFFSET(trialNum,channelLoop,mouseLoop)=tauValX;
%             
%             plot(f,x,exampleTrace)
%             
%             yL = get(gca,'YLim');
% 
%             line([tauValX tauValX],[yL(1) yL(2)],'Color','r');
%             line([tauValXCoeff tauValXCoeff],[yL(1) yL(2)],'Color','k');
%             
%             
%         end
%         
%         xlabelTitle=strcat(channelList{channel},'  OFFSET trial',num2str(trialNum));
%         
%         xlabel(xlabelTitle)
%         
%         
%         
%     end
%     
% end
% 
% %%
% 
% tauVectMatLabOFFSETReshape=[tauVectMatLabOFFSET(:,:,1);tauVectMatLabOFFSET(:,:,2);tauVectMatLabOFFSET(:,:,3)]/10;
% 
% %%
% 
% meantauVectMatLabOFFSET=mean(tauVectMatLabOFFSETReshape);
% 
% errTauVectMatLabOFFSET=squeeze(nanstd(tauVectMatLabOFFSETReshape,[],1)/sqrt(size(tauVectMatLabOFFSETReshape,1)));
% 
% %%
% 
% figure;
% 
% 
% for stimLoop=1:6;
% 
% stim=stimLoop;
% plot(stimLoop,tauVectMatLabOFFSETReshape(:,stim),'.m')
% hold on;
% errorbar(stimLoop+0.25,meantauVectMatLabOFFSET(:,stim),errTauVectMatLabOFFSET(:,stim),'om')
% 
% end
% 
% %%
% 
% stim=4;
% plot(2,tauVectMatLabOFFSETReshape(:,stim),'.g')
% hold on;
% errorbar(2.25,meantauVectMatLabOFFSET(:,stim),errTauVectMatLabOFFSET(:,stim),'og')
% 
% stim=5;
% plot(3,tauVectMatLabOFFSETReshape(:,stim),'.b')
% hold on;
% errorbar(3.25,meantauVectMatLabOFFSET(:,stim),errTauVectMatLabOFFSET(:,stim),'ob')
% 
% xlabel('offsets')
% ylabel('time(s)')
% 
% 
% saveas(gcf,'matlabTauOFFSETs','svg');
% saveas(gcf,'matlabTauOFFSETs','png');
% 
% %%  does matlab exponents exp1 or 2 for ONSET
% %
% 
% for trialNum=1:size(dFoverFBaselineOffsetMatDetrendFirst,2);
%     
%     for channelLoop=1:size(dFoverFBaselineOffsetMatDetrendFirst,3)
%         
%         figure('Renderer', 'painters', 'Position', [1+(channelLoop*280) 150 250 450]);
%         
%         channel=channelLoop;
%         
%         numMice=size(dFoverFBaselineOffsetMatDetrendFirst,4);
%         
%         for mouseLoop=1:numMice
%             
%             subplot(numMice,1,mouseLoop)
%             
%             %[val idx]=max(dFoverFBaselineOffsetMatDetrendFirst(:,trialNum,channel,mouseLoop));
%             
%             idx=offValVectON(channelLoop); % gives from offset of stim.
%             
%             exampleTrace=[];
%             
%             exampleTrace=dFoverFBaselineOffsetMatDetrendFirst(stimStart:idx,trialNum,channel,mouseLoop);
%             
%             x=1:1:size(exampleTrace,1);
%             
%             x=x';
%             
%             f=fit(x,exampleTrace,'exp1'); %exp1
%             
%             coeffVals=coeffvalues(f)
%             tauValXcoeff=coeffVals(2); % gets from matlab caluclated values.
%             
%             yVals = f(x);
%             tauValY=max((yVals)/100)*63.2;
%             
%             [~,tauValX] = (min(abs(yVals - tauValY)));  %tauValY
%             tauVectMatLabONSET(trialNum,channelLoop,mouseLoop)=tauValX;
%             
%             plot(f,x,exampleTrace)
%             
%             yL = get(gca,'YLim');
% 
%             line([tauValX tauValX],[yL(1) yL(2)],'Color','r');
%             line([tauValXcoeff tauValXcoeff],[yL(1) yL(2)],'Color','k');
%             
%         end
%         
%         xlabelTitle=strcat(channelList{channel},'  ONSET trial',num2str(trialNum));
%         
%         xlabel(xlabelTitle)
%         
%         
%     end
%     
% end
% 
% 
% %%
% 
% tauVectMatLabONSETReshape=[tauVectMatLabONSET(:,:,1);tauVectMatLabONSET(:,:,2);tauVectMatLabONSET(:,:,3)]/10;
% 
% %%
% 
% meantauVectMatLabONSET=mean(tauVectMatLabONSETReshape);
% 
% errTauVectMatLabONSET=squeeze(nanstd(tauVectMatLabONSETReshape,[],1)/sqrt(size(tauVectMatLabONSETReshape,1)));
% 
% %%
% 
% figure;
% 
% stim=3;
% plot(1,tauVectMatLabONSETReshape(:,stim),'.m')
% hold on;
% errorbar(1.25,meantauVectMatLabONSET(:,stim),errTauVectMatLabONSET(:,stim),'om')
% 
% stim=4;
% plot(2,tauVectMatLabONSETReshape(:,stim),'.g')
% hold on;
% errorbar(2.25,meantauVectMatLabONSET(:,stim),errTauVectMatLabONSET(:,stim),'og')
% 
% stim=5;
% plot(3,tauVectMatLabONSETReshape(:,stim),'.b')
% hold on;
% errorbar(3.25,meantauVectMatLabONSET(:,stim),errTauVectMatLabONSET(:,stim),'ob')
% 
% xlabel('onsets')
% ylabel('time(s)')
% 
% 
% saveas(gcf,'matlabTauONSETs','svg');
% saveas(gcf,'matlabTauONSETs','png');
% 

