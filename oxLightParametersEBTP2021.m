%% values to determine window clipping

stimOnset=300;
stimOffset=630;

%% values for filtering

linearFitType=1; % 0 takes from baseline, 1 takes from
order = 1;
framelen = 15; 

%% determine clipping, ranges to take values from, ranges to detrend from

stimStart=300;
stimStop=629;

baseStart=200;
baseStop=299;

stimValStart=stimStart; 
stimValStop=599;

detrendWindow=299;
baseDetrendStart1=1; %1  for detrending based on whole trace
baseDetrendStop1=detrendWindow; %299

baseMatVal=1; %0 = mean, 1= median, 2=max 3=STD.
stimMeanVal=0; %0 = mean, 1= median, 2=max.

%% decide which values you want to plot
 
getMeanOrMaxVals=1; %0 gives mean of stim window, 1 gives max of stim window.
 
meanOrMedianOfMice=0; %0 gives mean, 1 gives median.

%% for plotting

plotColList={'-k','-r','-m','-g','-b','-c'};

plotColListVal={'ok','or','om','og','ob','oc'};

forceAxis=1; %if 0, does not force axis, if 1, forces axis to yLForce



yLForceTraces=[-10 45]; %[-10 90];

yLForceVals=[-10 90]; %[-10 90];

