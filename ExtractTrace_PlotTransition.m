%% upload output of suite2p, plot ROI traces, FOV z-projection, ROIs coordinates and then save traces in  one .m file per FOV.

close all
clear all

%%%%%%% choose FOV to analyse, based on the list of FOVs in the FOVs_used function
oxlight=1; %set to 0 if it is an oxlight control, set to 1 if it is oxlight1 mouse.
FOV_toanalyse=7;
[img_root, img_file_a, img_file_tosave, mouseID]=FOVs_used(oxlight,FOV_toanalyse);
%%%%%%%

asleeptoawake=1; %1 if mouse starts asleep and ends awake, 0 if mouse starts awake and ends asleep
baseline_type=1; % 0 if you don't want to calculate df/f, 1 for df/f0, with f0 being the median of the asleep period
one_roi=1; % switch between 1 to analyse the whole FOV, or 0  for individual ROIs.

if oxlight
    savepath=(['C:\Users\oxlight\gradual_transitions']);
else
    savepath=(['C:\Users\ctrl\gradual_transitions']);
end

load(['C:\Users\suite2p\plane0\Fall.mat'])
frate=ops.fs; %imaging frame rate

%df/f traces (or raw traces if baseline_type==0
ROIs_corrected = [];
dFoF = [];

%only analyse fluorescence from pixels clusters identified as neurons
list_ROIs=find(iscell(:,1)==1); %indeces of the pixels clusters identified as cells
F = F(list_ROIs,:);

traces_plot=figure(10);clf
colors=colormap(parula(size(F,1)));
for rr=1:size(F,1)
    
    ROIs_corrected = [ROIs_corrected; (F(rr,:))];
    if baseline_type==0
        dFoF = [dFoF; ROIs_corrected(rr,:)];
    elseif baseline_type==1
        [transition, dFoF_temp,F0, F0_std]=baseline_bistate(ROIs_corrected(rr,:),frate,asleeptoawake);
        dFoF = [dFoF; dFoF_temp];
    end
    
    hold on
    if one_roi
        plot(1:length(dFoF),dFoF(rr,:), 'LineWidth',0.5, 'Color', colors(rr,:))
    else
        plot(1:length(dFoF),dFoF(rr,:)+rr*5, 'LineWidth',0.5, 'Color', colors(rr,:)) %use this when plotting more than one ROIs
    end
end

ylabel('df/f')
xlabel('imaging frames')
% xlim([0 length(F)])

transitions_tosave=[savepath '/' mouseID ' _' img_root ' _' img_file_tosave '_gradual-transitions.mat'];
save(transitions_tosave, 'F0', 'transition', 'frate', 'dFoF')


%calculate and then save mean dfof and std for each trace
avg_dFoF=mean(dFoF');
std_dFoF=std(dFoF');

%heatmap traces
figure(14);clf
imagesc(dFoF)
ylabel('ROI')
xlabel('imaging frames')

% identify neurons' coordinates
if one_roi==0
    shape_rois=zeros(size(ops.meanImg,1), size(ops.meanImg,2));
    stat=stat(list_ROIs);
    for gg=1:length(stat)
        coords(gg,:)=stat{gg}.med; %median x and y coordinates for each good neuron in the field. First column is y, second column is x!!
        for xx= 1:length(stat{gg}.xpix)
            shape_rois(stat{gg}.ypix(xx), stat{gg}.xpix(xx))= 1;
        end
    end
    
    figure (12);clf
    imagesc(shape_rois)
    title([num2str(mouseID) ' ' num2str(img_root) ' '  num2str(img_file_a)])
    for ww=1:length(coords)
        text(coords(ww,2), coords(ww,1),num2str(ww) ,'Color',[1 0 0])
    end
end

%% analysis on traces including the whole FOV. Plot each FOV individually, and the mean trace across FOVs. Do so for both oxlight FOVs and ctrl FOVs.

clearvars -except asleeptoawake
close all

%the roots here are the locations of the  files saved while running the
%"plot_traces" script. OxLight1 and OxLight-ctr are in two separate folders.
root_oxlight="C:\Users\oxlight\gradual_transitions\";
root_ctrl="C:\Users\ctrl\gradual_transitions\";

if asleeptoawake==1 %analyse asleep to awake
    transistion_files_oxlight = dir(fullfile(root_oxlight + "*_asleeptoawake*.mat")); %list all files saved while using the "plot_traces" script
elseif asleeptoawake==0 %analyse reverse
    transistion_files_oxlight = dir(fullfile(root_oxlight + "*_reverse_*.mat"));
end
allmice_transition_oxlight=[];
allmice_downsampled_oxlight=[]; %downsample traces just for plotting purposes. The reason is because FOVs were acquired with different frame rates.

all_durations=[];
all_frate=[];
for xx = 1:length(transistion_files_oxlight)
    data = load(transistion_files_oxlight(xx).name);
    
    all_durations=[all_durations max((linspace(1,length(data.dFoF),length(data.dFoF)))./data.frate*1)];
    all_frate=[all_frate data.frate];
end
clear data

colorstouse=colormap(jet(length(transistion_files_oxlight)));
for ii = 1:length(transistion_files_oxlight)
    sampled_points=[];
    downsampled_points=[];
    data = load(transistion_files_oxlight(ii).name);
    transistion_files_oxlight(ii).name
    data.dFoF_cut=data.dFoF(1:(length(data.dFoF)*min(all_durations))/all_durations(ii));
    sampled_points=linspace(1,length(data.dFoF_cut),length(data.dFoF_cut))./data.frate*1;
    downsampled_points=linspace(1/min(all_frate),max(sampled_points),round(max(sampled_points)*min(all_frate)));
    data.dFoF_downsampled=interp1(sampled_points, data.dFoF_cut, downsampled_points);
    length(data.dFoF_downsampled)
    
    figure(1);hold on
    plot((linspace(1,length(data.dFoF),length(data.dFoF)))./data.frate*1, data.dFoF, 'Linewidth', 0.5, 'Color', colorstouse(ii, :))
    
    figure(10); hold on
    plot((linspace(1,length(data.dFoF_downsampled),length(data.dFoF_downsampled)))./min(all_frate)*1, data.dFoF_downsampled, 'Linewidth', 0.4,'Color',[0 1 0])% colorstouse(ii, :))
    
    allmice_transition_oxlight=[allmice_transition_oxlight data.transition];
    if ii == 1 %run this in order to add the average trace across FOVs to the plot
        set_traces_length = length(data.dFoF_downsampled);
    else
        set_traces_length = size(allmice_downsampled_oxlight,2);
    end
    allmice_downsampled_oxlight=[allmice_downsampled_oxlight; data.dFoF_downsampled(1:set_traces_length)];
    
end
mean_acrossoxlight_downsampled=mean(allmice_downsampled_oxlight);

figure(1)
ylim([-0.2 0.7])
xlabel('time (s)')
ylabel ('df/f - oxlight1')
set(gca,'TickDir', 'out', 'FontSize', 14, 'Linewidth', 2.5)
box off
figure(10); hold on
plot((linspace(1,length(mean_acrossoxlight_downsampled),length(mean_acrossoxlight_downsampled)))./min(all_frate)*1, mean_acrossoxlight_downsampled, 'g', 'Linewidth', 1.5)
ylim([-0.2 0.9])
xlim([0 360])

xlabel('time (s)')
ylabel ('df/f - oxlight1')
set(gca,'TickDir', 'out', 'FontSize', 14, 'Linewidth', 2.5)

%run the same operations, but for oxlight-ctr FOVs
clearvars -except allmice_transition_oxlight root_ctrl asleeptoawake

if asleeptoawake==1
    transistion_files_ctrl = dir(fullfile(root_ctrl + "*_asleeptoawake*.mat"));
elseif asleeptoawake==0
    transistion_files_ctrl = dir(fullfile(root_ctrl + "*_reverse_*.mat"));
end

allmice_transition_ctrl=[];
allmice_downsampled_ctrl=[];

all_durations_ctrl=[];
all_frate_ctrl=[];
for xx = 1:length(transistion_files_ctrl)
    data = load(transistion_files_ctrl(xx).name);
    
    all_durations_ctrl=[all_durations_ctrl max((linspace(1,length(data.dFoF),length(data.dFoF)))./data.frate*1)];
    all_frate_ctrl=[all_frate_ctrl data.frate];
end
clear data

colorstouse_ct=colormap(jet(length(transistion_files_ctrl)));
for i = 1:length(transistion_files_ctrl)
    sampled_points_ctrl=[];
    downsampled_points_ctrl=[];
    data = load(transistion_files_ctrl(i).name);
    transistion_files_ctrl(i).name
    
    data.dFoF_cut_ctrl=data.dFoF(1:(length(data.dFoF)*min(all_durations_ctrl))/all_durations_ctrl(i));
    sampled_points_ctrl=linspace(1,length(data.dFoF_cut_ctrl),length(data.dFoF_cut_ctrl))./data.frate*1;
    downsampled_points_ctrl=linspace(1/min(all_frate_ctrl),max(sampled_points_ctrl),round(max(sampled_points_ctrl)*min(all_frate_ctrl)));
    data.dFoF_downsampled_ctrl=interp1(sampled_points_ctrl, data.dFoF_cut_ctrl, downsampled_points_ctrl);
    length(data.dFoF_downsampled_ctrl)
    
    figure(2);hold on
    plot((linspace(1,length(data.dFoF),length(data.dFoF)))./data.frate*1, data.dFoF, 'Linewidth', 0.5, 'Color', colorstouse_ct(i, :))
    
    figure(20); hold on
    plot((linspace(1,length(data.dFoF_downsampled_ctrl),length(data.dFoF_downsampled_ctrl)))./min(all_frate_ctrl)*1, data.dFoF_downsampled_ctrl, 'Linewidth', 0.4,'Color', [0 0 0])% colorstouse_ct(i, :))
    
    allmice_transition_ctrl=[allmice_transition_ctrl data.transition];
    if i == 1
        set_traces_length = length(data.dFoF_downsampled_ctrl);
    else
        set_traces_length = size(allmice_downsampled_ctrl,2);
    end
    allmice_downsampled_ctrl=[allmice_downsampled_ctrl; data.dFoF_downsampled_ctrl(1:set_traces_length)];
    
end
mean_acrossctrl_downsampled=mean(allmice_downsampled_ctrl);

figure(2)
ylim([-0.2 0.7])
xlabel('time (s)')
ylabel ('df/f - oxlight ctr')
set(gca,'TickDir', 'out', 'FontSize', 14, 'Linewidth', 2.5)
figure(20); hold on
plot((linspace(1,length(mean_acrossctrl_downsampled),length(mean_acrossctrl_downsampled)))./min(all_frate_ctrl)*1, mean_acrossctrl_downsampled, 'k', 'Linewidth', 1.5)
ylim([-0.2 0.7])
xlim([0 420])
xlabel('time (s)')
ylabel ('df/f - oxlight-ctr')
set(gca,'TickDir', 'out', 'FontSize', 14, 'Linewidth', 2.5)


avg_oxlight=mean(allmice_transition_oxlight);
avg_ctrl=mean(allmice_transition_ctrl);

%% plot values of df/f change from 1st to last minute of imaging, in both ctrls and oxlight FOVs. Plot both mean across FOVs and individual FOVs.

figure(3)
errorbar(1,avg_oxlight, std(allmice_transition_oxlight)./sqrt(length(allmice_transition_oxlight)),'Linewidth', 2.5, 'Color', [0 0 0])
hold on
errorbar(1.3,avg_ctrl, std(allmice_transition_ctrl)./sqrt(length(allmice_transition_ctrl)),'Linewidth', 2.5, 'Color', [0 0 0])
hist_data = [avg_oxlight avg_ctrl];
bar(hist_data, 'k')
hold on
scatter((ones(1,length(allmice_transition_oxlight))).*1.25, [allmice_transition_oxlight], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [.7 .7 .7])
hold on
scatter((ones(1,length(allmice_transition_ctrl))).*2.25, [allmice_transition_ctrl], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [.7 .7 .7])
ylim([-0.2 0.6])
xlim([0.9 1.35])

box off
set(gca,'TickDir', 'out', 'Linewidth', 2, 'xtick', [1 2], 'xlim', [0.5 2.5], 'XTickLabel', {'oxLight', 'control'}, 'FontSize', 14)
ylabel('df/f')

% stats comparing ctrl and oxlight
[h,p] = ttest2(allmice_transition_oxlight, allmice_transition_ctrl)
[p,h] = ranksum(allmice_transition_oxlight, allmice_transition_ctrl)
