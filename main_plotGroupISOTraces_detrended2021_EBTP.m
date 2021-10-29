%% Dose response FREQUENCY

%update should contain plotOxLightDoseResponseDETRENDfixed

clear

%% load dose response matrix

oxOrCtrl=1;  %0 loads oxlight, 1 loads ctrl, 2 loads 10Hz train data

%%

if oxOrCtrl==0; 

fileLoadName=dir('mainIsoDoseResponseMatOX.mat');  %mainIsoDoseReponseMat10HzTrain.mat

titleStub='oxLight';

elseif oxOrCtrl==1; 

fileLoadName=dir('mainIsoDoseResponseMatCTRL.mat');  %mainIsoDoseReponseMat10HzTrain.mat    

titleStub='oxLightCtrl';

elseif oxOrCtrl==2; 

fileLoadName=dir('mainIsoDoseResponseMat10HzTrain.mat');  %mainIsoDoseReponseMat10HzTrain.mat    
    
titleStub='oxLight10Hz';

end
    
fileLoadName=fileLoadName.name;
load(fileLoadName);

%% input Values

fiveHzChannel=1; %1 for control mice, %3 for normal oxight

oxLightParametersEBTP2021

legendList={'catchTrial','1Hz','5Hz','10Hz','20Hz','las onset','las offset'};

%% clip response Matrix to size

midPoint=ceil(size(mainIsoDoseReponseMat,1)/2);

mainIsoDoseReponseMat=mainIsoDoseReponseMat(midPoint-stimOnset:midPoint+((stimOffset-stimOnset)*3),:,:,:);

%% filter traces

lengthMast=size(mainIsoDoseReponseMat,1);

    
    for mouseLoop=1:size(mainIsoDoseReponseMat,4)
        
        tempMouseFiltMat=squeeze(mainIsoDoseReponseMat(:,:,:,mouseLoop));
        
        for trialLoop = 1:size(tempMouseFiltMat,2);
            
            clear sgTempFiltMat
            
            tempTrialFiltMat=squeeze(tempMouseFiltMat(:,trialLoop,:));
            
            sgTempFiltMat = sgolayfilt(tempTrialFiltMat,order,framelen);
            
            %sgTempFiltMat = reshape(sgTempFiltMat,3,5)
            
            trialFiltMat(:,trialLoop,:)=sgTempFiltMat;
            
        end
        
        filtMouseMat(:,:,:,mouseLoop)=trialFiltMat;
        
    end
    
    
    %% plot sample of filtered and unfiltered traces
    
    figure; %plot sample of filtered and unfiltered traces
    
    numSubplots=size(mainIsoDoseReponseMat,3);
    
    for subPlotLoop=1:numSubplots
        
        subplot(numSubplots,1,subPlotLoop)
        
        plot(mainIsoDoseReponseMat(:,1,subPlotLoop,1),'-k','lineWidth',2)
        hold on
        plot(filtMouseMat(:,1,subPlotLoop,1),'-g')
        
        
    end
    
    legend('raw','sgolay')
    
    mainIsoDoseReponseMat=filtMouseMat;
    

%% get means of trials per animal, then means across animal.

meanRawMainTrial=squeeze(mean(mainIsoDoseReponseMat,2));
meanRawMainTrialProb=squeeze(mean(meanRawMainTrial,3));

%% detrend and plots each trace

tempdFoverFMatrix=mainIsoDoseReponseMat;

plotOxLightDoseResponseDETRENDfixed

rawDETRENDED=mainTrialDetrendMat;

save('rawDETRENDED', 'rawDETRENDED' );

%%

plotOxLightDoseResponseDetrendFirstThenZscore2021TPEB

