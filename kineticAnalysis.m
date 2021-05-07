%%
clear all;
[tau,latency] = sensorKinetics();
%input arguments:
%   name: name of mat file generated by rawdataK_ON.py or rawdataK_OFF.py
%   interval: frametime
%   frameT0: beginning of the experiment, manually determined by
%   experimenter
%   frameTend: end of the experiment, manually determined by experimenter
%   ON/OFF: Experiment type (Agonist = 0, Antagonist = 1)
%   temporal resolution: =1 if further smoothing of the data is needed


function [Tau,Latency] = sensorKinetics(name,interval,FrameT0,FrameTend,Exp,Avg)

experimentName = name; %Experiment name
frameTime = interval; %Time per frame
aa = FrameT0;%% Start
bb = FrameTend; %% End

load(strcat(experimentName,'.mat'));
fs = round(1000/frameTime);
fpass = 1; 

%%read data
TexRedPix1 = double(TexRedPix(:,aa:bb));
GFP1 = double(GFP(:,aa:bb));

[Nrow,Ncolumn] = size(GFP1);
GdFF = zeros(Nrow,Ncolumn);
RdFF = zeros(Nrow,Ncolumn);

%%denoising
for i = 1:Nrow
    GdFF(i,:) = lowpass(GFP1(i,:)/mean(GFP1(i,1:50)),fpass,fs);
end
Gmax = mean(GdFF(:,end-100:end-50),2);
GnormValue = GdFF(:,1:end-50)./Gmax;

for i = 1:Nrow
    RdFF(i,:) = lowpass(TexRedPix1(i,:)/mean(TexRedPix1(i,1:50)),fpass,fs);
end
Rmax = mean(mean(RdFF(:,end-100:end-50)));
RnormValue = RdFF(:,1:end-50)/Rmax;   

X = (0:1:size(GnormValue,2)-1);

plot(X*frameTime,GnormValue,'color',[0.85,0.85,0.85]);
hold on;


%%TexRed signal latency
DR1R = gradient(mean(RnormValue));
timeZero = find(DR1R == max(DR1R));
Rplateau = mean(RnormValue(:,end-50:end));
MEANR = mean(Rplateau);
RBOOLEAN = (mean(RnormValue) < MEANR*0.85);
RRRR =tabulate(RBOOLEAN);
Latency = (RRRR{2,2} - timeZero)*frameTime;

%adjusting temporal resolution
if Avg == 0
    plot(X*frameTime,GnormValue,'color',[0.85,0.85,0.85]);
    hold on;
else
    GnormValue = GnormValue(:,1:end-1);
    iNterval = 10;
    Nint = size(GnormValue,2)/iNterval; 
    Gnorm1 = reshape(GnormValue,Nrow,iNterval,Nint);
    GnormValue = squeeze(mean(Gnorm1,2));
    X = (0:1:size(GnormValue,2)-1).*iNterval;
    timeZero = round(timeZero./iNterval);
    plot(X*frameTime,GnormValue,'color',[0.85,0.85,0.85]);
    hold on;
end

%fitting
if Exp == 0
    [Tau,Yplot] = onePhaseAssoFitting(timeZero,X*frameTime,mean(GnormValue));
    hold on;
    plot(X*frameTime+timeZero*frameTime*iNterval,Yplot,'LineWidth',2);
else
    [Tau,Yplot] = onePhaseDissoFitting(timeZero,X*frameTime,mean(GnormValue));
    hold on;
    plot(X*frameTime+timeZero*frameTime,Yplot,'LineWidth',2);
    Latency = 0;
end
end

function [TauOn,Yvalue] = onePhaseAssoFitting(T0,XData,DataSet)


Y0 = mean(DataSet(1:T0));
Ymax = mean(DataSet(end-1:end));
KK = 1e-4;

onePhaseAsso = @(k,x)((Ymax-Y0).*(1-exp((-k).*(x))));

beta = nlinfit(XData,DataSet,onePhaseAsso,KK);
Yvalue = onePhaseAsso(beta,XData)+Y0;
TauOn = 1/beta;

end

function [TauOFF,Yvalue] = onePhaseDissoFitting(T0,XData,DataSet)


Y0 = mean(DataSet(1:T0));
Ymax = mean(DataSet(end-1:end));
KK = 1e-4;
onePhaseDisso = @(k,x)(Ymax+(Y0-Ymax).*(exp((-k).*(x))));

beta = nlinfit(XData,DataSet,onePhaseDisso,KK);
Yvalue = onePhaseDisso(beta,XData);
TauOFF = 1/beta;

end