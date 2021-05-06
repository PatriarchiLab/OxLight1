function [conf] = DFFanalysisConfig()

%% Input/output options %%
conf.DFFanalysisPath = 'Y:\Luca\MatlabCode\DFFanalysis\_DFFanalysis'; % Path to the main FRETanalysis folder (the one containing the core code)
conf.inputFilePath = ''; % Full path to the input file
conf.inputFileType = 'Tiff'; % Options: Tiff, ZeissZen
conf.outName = ''; % Name of the current analysis (i.e. name of the subdir in Results)

%% Analysis options %%
%% IMPORTANT: All operations are sequential, and they happen in the same order as they are listed below. %%

% Channel selection
conf.chInfo.totChannelsInInput = 1; % Total number of channels in input image (including those not relevant for FRET)
conf.chInfo.channelToUse = 1; % Select the channel to use for DF/F 

% Spatial resolution reduction or smoothing
conf.spaceBinning = 0; % SWITCH (0/1): reduce the image resolution by a factor conf.downsizeFactor
conf.spaceBinningFactor = 1; % Reduce the resolution by this factor (e.g. 512x512 with factor 2 becomes 256x256, with factor 4 becaomes 128x128); 
conf.gaussSmooth = 0; % SWITCH (0/1): spatially smooth the images (gaussian smoothing)
conf.gaussSmoothStdDev = 1; % Standard deviation of the gaussian used for smoothing

% Tresholding 
conf.minThresh = 0; % SWITCH (0/1): apply minimun value threshold
conf.minThreshValue = 300; % Minimum threshold value
conf.maxThresh = 0; % SWITCH (0/1): apply maximum value threshold
conf.maxThreshValue = 2000; % Maximum threshold value
conf.ADVthresh = 0; % SWITCH (0/1): use an advanced threshold method 
conf.ADVthreshType = 'Li'; % Type of advanced threshold. Options: 'Otsu', 'Li', 'Chastagnier'

% Ratio type and thresholding
conf.unboundFrames = XX:XX; % Which frames of the image correspond to the unbound sensor
conf.boundFrames = XX:XX; % Which frames of the image correspond to the bound sensor
conf.ratioThresholding = 0; % SWITCH (0/1): turn into nans all ratio values outside a defined interval
conf.ratioThresholdingLimits = [0 15]; % The defined interval ([min max]) for histogram-based denoising

%% Output images options %%
conf.outImages.print = 1; % SWITCH (0/1): print tiff images
conf.outImages.printMask = 1; % SWITCH (0/1): print the segmentation/tresholding mask, as binary image (depends on conf.outImages.combine)
conf.outImages.colormap = 'kvroyw'; % See colorMapsDoc file for list of available colormaps
conf.outImages.cmapMin = 0; % Minimum value for the colormap (use -1 to autodetect minimum)
conf.outImages.cmapMax = 15; % Maximum value for the colormap (use -1 to autodetect maximum)
conf.outImages.highlightSaturation = 0; % Color white all pixels above Imax
conf.outImages.printCB = 1; % SWITCH (0/1): print the colorbar as separate image
conf.outImages.printSteps = 1; % SWITCH (0/1): print a Matlab Figure with all steps of intensity weighting

%% ANYTHING BELOW THIS LINE IS FOR INTERNAL WORKING OF THE CODE. DO NOT CHANGE!            
conf.configPath = strcat(mfilename('fullpath'), '.m');
                                         
end

