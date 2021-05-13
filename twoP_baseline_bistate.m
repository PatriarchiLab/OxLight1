function [transition, dF_bas, F0, F0_std]=twoP_baseline_bistate(input,frate,asleeptoawake)

% m.panniello 2020
% calculate (f-f0)/f0. F0 is the median of the pixel distribution during the anaesthetized period.


secs_long=length(input)/frate; %how long is the acquisition in seconds

one_minute=floor((length(input)*60)/secs_long); %how many frames correspond to one minute of the t-series.

if asleeptoawake==1
    F0=median(input(1:one_minute));
    F0_std=std(input(1:one_minute));
else
    F0=median(input(one_minute*2:end));
    F0_std=std(input(one_minute*2:end));
end

dF_bas=(input-F0)/F0;

if asleeptoawake==1
    startoftrace_fluo=median(dF_bas(1:floor(one_minute)));
    endoftrace_fluo=median(dF_bas(end-floor(one_minute):end));
else
    startoftrace_fluo=median(dF_bas(1:floor(one_minute)*2));
    endoftrace_fluo=median(dF_bas(end-floor(one_minute)*2:end));
end


transition=endoftrace_fluo-startoftrace_fluo; %quantify progressive increase or decrease of signal following the asleep to awake (or awake to asleep) transition of the mouse
