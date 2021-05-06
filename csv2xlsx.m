%% N.B. 2019/11/06
%For .csv file exported from fiji
%Exports two .xlsx files: 
%   *_dF.xlsx for calculated dF data
%   *_average.xlsx for maximum response dF
    
%frameNo 1 - n: index of event frame (e.g. addition of drug)


%% load data directory
clear all;


dataTank = dir('*.csv'); 
n = size(dataTank,1);
frameNo1 = []; %REQUIRES MANUAL IMPUT:


%% importing csv data

for i = 1:n
    
    cd('C:\Users\lduffe\Desktop\ICL2_Mutants analysis')
    fileName = dataTank(i).name; 
    rawData = table2array(readtable(fileName));    
    diffBase = frameNo1(i) - min(frameNo1)+1;
    
    
    for j = 2:size(rawData,2)
        alignedData(:,j-1) = rawData(diffBase:end,j);
        
        %getting dF
        %column one = frame number i.e. time
        dF(:,1) = [1:size(alignedData,1)];
        dF(:,j) = alignedData(:,j-1) / mean(alignedData(1:min(frameNo1),j-1)) - 1;                    
    end
    
    %% averaging last ten frames - for current experiment only
     l = size(dF,1);
     w = size(dF,2);
        for k = 2:w
            averagedF(k-1,1) = mean(dF(l-12:l,k));
            transaveragedF = averagedF'
        end
        
    %% save current workspace n data
    mkdir('C:\Users\lduffe\Desktop\ICL2_Mutants analysis\Processed_data_temp')
    cd('C:\Users\lduffe\Desktop\ICL2_Mutants analysis\Processed_data_temp')
    save(erase(fileName,'.csv'),'rawData','truncData','dF','averagedF','transaveragedF');
    writetable(array2table(dF), strcat(erase(fileName,'.csv'),'_dF','.xlsx'));
    %this gives you the dF data in an xlsx filetype
    %% writetable(array2table(averagedF), strcat(erase(fileName,'.csv'),'_average','.xlsx'));
    %this gives you the averaged dF of each ROI, ignore the first column in
    %the xlsx file
    clear rawData;
    clear truncData;
    clear averagedF
    clear dF;
end

        
        
        
        
        
        