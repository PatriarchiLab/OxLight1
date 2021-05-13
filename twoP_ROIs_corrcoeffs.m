% use this script to:
% 1)identify, in each FOV separately, the one minute of highest activity (most active minute in Duffet et al)
% 2)plot traces from ACTIVE ROIs and from INACTIVE ROIs (both raw activity and deconvolved activity)
% 3)calculate pairwise pearson's correlation coefficient on each ACTIVE ROI pair and on each INACTIVE ROI pair
% 4)save average corr coeff across pairs for each FOV, then plot these for all FOVs, in oxlight and ctrl FOVs.

%% run this section for each individual FOV
close all
clear all

%%%%%%% choose FOV to analyse, based on the list of FOVs in the FOVs_used function
oxlight=1; %set to 0 if it is an oxlight control, set to 1 if it is oxlight mous
FOV_toanalyse= 7;
[img_root, img_file_a, img_file_tosave, mouseID]=twoP_FOVs_used(oxlight,FOV_toanalyse);

% set up the following variables:
after_InactiveRois = 11; %rois with high deconvolution value (active rois) start from this ROI index in each FOV;
plot_paper=1; %plot data for manuscript

if oxlight
    savepath=(['\home\results\oxlight\active_rois']);
    savepath_inactive=(['\home\results\oxlight\inactive_rois']);
    root_deconv=(['\home\results\oxlight\deconvolution']);
    root_fullFOVtrace=(['\home\results\oxlight\gradual_transitions']);
    
else
    savepath=(['\home\results\ctrl\active_rois']);
    savepath_inactive=(['\home\results\ctrl\inactive_rois']);
    root_deconv=(['\home\results\ctrl\deconvolution']);
    root_fullFOVtrace=(['\home\results\ctrl\gradual_transitions']);
end

load([root_deconv '/' mouseID ' _' img_root ' _' img_file_tosave '_deconvolved.mat'])
load([root_fullFOVtrace '/' mouseID ' _' img_root ' _' img_file_tosave '_gradual-transitions.mat'])
fullFOVtrace=dFoF; %this is the dFoF trace calculated in the ExtractTrace_PlotTransition script, using one ROI including the whole FOV.
clear dFoF F0
load(['C:\Users\suite2p\plane0\Fall.mat'])

%how many frames correspond to one minute of the series.
secs_long=length(fullFOVtrace)/frate;
one_minute=floor((length(fullFOVtrace)*60)/secs_long);

%apply moving window to find one minute showing highest fluorescent variations (referred to as "most active minute" in Duffet et al)
search_most_active=one_minute*2:one_minute*6; %only look for the most active minute between the 3rd and 7th minute of imaging
fullFOVtrace_awakening=fullFOVtrace(search_most_active);
[max_std, index_totake]=max(movstd(fullFOVtrace_awakening,one_minute));
most_active_minute=fullFOVtrace_awakening(index_totake:index_totake+one_minute);
figure(100)
plot(most_active_minute)

%separately analyse fluorescence from pixels clusters identified as active ROIs and inactive ROIs
dFoF = [];
dFoF_inactive = [];
list_ROIs=find(iscell(:,1)==1); %indeces of the pixels clusters identified as ROIs
F = F(list_ROIs,:);

fig_traces_plot=figure(1);clf
if oxlight
    colors=[0 0.63 0];
else
    colors=[0 0 0];
end

%plot one or more ACTIVE rois
for rr=after_InactiveRois:size(F,1)
    
    dFoF = [dFoF; F(rr,:)];
    
    hold on
    plot(1:length(dFoF),dFoF(rr-10,:)-rr*1, 'LineWidth',0.5, 'Color', colors)
end
ylabel('active ROI index')
xlabel('imaging frames')
title('raw signal active rois')
hold on
%add dotted lines indicating the most active minute
plot([one_minute*2+index_totake one_minute*2+index_totake],[min(dFoF,[],'all') max(dFoF,[],'all')],'--k')
hold on
plot([one_minute*2+index_totake+one_minute one_minute*2+index_totake+one_minute],[min(dFoF,[],'all') max(dFoF,[],'all')],'--k')
xlim([1 length(dFoF)])


%plot one or more INACTIVE rois
figure(10)
for rr=1:(after_InactiveRois-1)
    
    dFoF_inactive = [dFoF_inactive; F(rr,:)];
    
    hold on
    plot(1:length(dFoF_inactive),dFoF_inactive(rr,:)-rr*1, 'LineWidth',0.5, 'Color', colors)
end
ylabel('inactive ROI index')
xlabel('imaging frames')
title('raw signal inactive rois')
hold on
%add dotted lines indicating the most active minute
plot([one_minute*2+index_totake one_minute*2+index_totake],[min(dFoF_inactive,[],'all') max(dFoF_inactive,[],'all')],'--k')
hold on
plot([one_minute*2+index_totake+one_minute one_minute*2+index_totake+one_minute],[min(dFoF_inactive,[],'all') max(dFoF_inactive,[],'all')],'--k')
xlim([1 length(dFoF_inactive)])

%heatmap traces
fig_heat_all=figure(2);clf; %set(gcf, 'Position', [629   447   347   314])
allrois_matrix=[dFoF_inactive; dFoF];
imagesc(allrois_matrix)
colorbar
colormap bone
ylabel('ROI')
xlabel('imaging frames')
title('heatmap df/f inactive and active rois')
hold on
plot([one_minute*2+index_totake one_minute*2+index_totake],[0 size(F,1)+1],'-k')
hold on
plot([one_minute*2+index_totake+one_minute one_minute*2+index_totake+one_minute],[0 size(F,1)+1],'-k')


% identify ROIs' coordinates (active and inactive) and plot their position (their centers) on the projection of the deconvolved FOV.
shape_rois=zeros(size(ops.meanImg,1), size(ops.meanImg,2));
stat=stat(list_ROIs);
for gg=1:length(stat)
    coords(gg,:)=stat{gg}.med; %median x and y coordinates for each roi in the field. First column is y, second column is x!!
    for xx= 1:length(stat{gg}.xpix)
        shape_rois(stat{gg}.ypix(xx), stat{gg}.xpix(xx))= 1;
    end
end
fig_coords=figure (3);clf; %set(gcf, 'Position', [82   257   547   387])
imagesc(avg_luc1)
colorbar
caxis([0 0.56])
title([num2str(mouseID) ' ' num2str(img_root) ' '  num2str(img_file_a)])
for ww=1:length(coords)
    if ww<11
        text(coords(ww,2), coords(ww,1),num2str(ww) ,'Color',[0 0 0],'FontSize',11)
    else
        text(coords(ww,2), coords(ww,1),num2str(ww) ,'Color',[1 0 0],'FontSize',11)
        
    end
end

%now identify rois on deconvolved t-series and plot the deconvolved activity (extended data figure 8, panel g)
fig_deconvtraces=figure(4); hold on; %set(gcf, 'Position', [978    60   835   318])
allrois_deconv_traces=[];
for rr=1:length(stat) %for each roi
    trace_allframes=[];
    for ff=1:size(luc1,3) %for each frame
        thisroisarea=[];
        for xx= 1:length(stat{rr}.xpix)
            this_roispxls = luc1(stat{rr}.ypix(xx), stat{rr}.xpix(xx),ff);
            thisroisarea=[this_roispxls thisroisarea];
        end
        mean_perframe=mean(thisroisarea);
        trace_allframes=[trace_allframes mean_perframe];
    end
    allrois_deconv_traces=[allrois_deconv_traces; trace_allframes]; %obtained deconvolved traces across the frames used for creating luc1
    plot(1:size(allrois_deconv_traces,2),allrois_deconv_traces(rr,:)-rr*1, 'LineWidth',0.5, 'Color', colors)
end
ylabel('ROI index')
xlabel('imaging frames')
title('deconvolved activity')
xlim([1 size(luc1,3)])
hold on
plot([index_totake index_totake],[0 -size(allrois_deconv_traces,1)],'-k')
hold on
plot([index_totake+one_minute index_totake+one_minute],[0 -size(allrois_deconv_traces,1)],'-k')

fig_heatdeconv=figure(5); %set(gcf, 'Position', [629    63   349   313])
imagesc(allrois_deconv_traces)
colorbar
title('heatmap deconvolved rois')
hold on
plot([index_totake index_totake],[0 size(allrois_deconv_traces,1)+1],'-k')
hold on
plot([index_totake+one_minute index_totake+one_minute],[0 size(allrois_deconv_traces,1)+1],'-k')

%%%%%% calculate correlation coefficients
cmbntns=nchoosek(1:after_InactiveRois-1,2);
%corr coeff between raw activity in pairs of ACTIVE rois, for the most active minute
allrois_corrcoef=[];
allrois_traces =dFoF(:,one_minute*2:one_minute*6);
for cc=1:length(cmbntns)
    data_1 = allrois_traces(cmbntns(cc,1),index_totake:index_totake+one_minute);
    data_2 = allrois_traces(cmbntns(cc,2),index_totake:index_totake+one_minute);
    data_1 = detrend(data_1,2,round(length(data_1)*[0.25 0.5 0.75 1]));
    data_2 = detrend(data_2,2,round(length(data_2)*[0.25 0.5 0.75 1]));
    
    r_coeff=corrcoef(data_1-mean(data_1),data_2-mean(data_2));
    
    allrois_corrcoef=[allrois_corrcoef r_coeff];
end
allrois_corrcoef=(allrois_corrcoef(2,:))';
allcoeffs_raw=allrois_corrcoef(allrois_corrcoef~=1); %average correlation coefficient (r, lag0) across all pairs, without autocorrelation
avg_corrcoef_raw=mean(allcoeffs_raw);

%plot correlation between ACTIVE ROI pairs
figure(15)
raw_traces_pearson = (allrois_traces(1:end,index_totake:index_totake+one_minute))';
raw_detrended = detrend(raw_traces_pearson,2,round(length(raw_traces_pearson)*[0.25 0.5 0.75 1]));
[R,PValue] = corrplot(raw_detrended);%,'tail','right');
fig_corrcoeff_raw=figure(16);
colormap(redblue)
imagesc(R)
caxis([-1 1])
colorbar
ylabel ('roi number')
xlabel ('roi number')
title('corr coeff on ACTIVE roi pairs')

clear data_1 data_2

%corr coeff between raw activity in pairs of INACTIVE rois, for the most active minute
allrois_corrcoef_inactive=[];
allrois_traces_inactive =dFoF_inactive(:,one_minute*2:one_minute*6);
for cc=1:length(cmbntns)
    data_1 = allrois_traces_inactive(cmbntns(cc,1),index_totake:index_totake+one_minute);
    data_2 = allrois_traces_inactive(cmbntns(cc,2),index_totake:index_totake+one_minute);
    data_1 = detrend(data_1,2,round(length(data_1)*[0.25 0.5 0.75 1]));
    data_2 = detrend(data_2,2,round(length(data_2)*[0.25 0.5 0.75 1]));
    
    r_coeff_inactive=corrcoef(data_1-mean(data_1),data_2-mean(data_2));
    
    allrois_corrcoef_inactive=[allrois_corrcoef_inactive r_coeff_inactive];
end
allrois_corrcoef_inactive=(allrois_corrcoef_inactive(2,:))';
allcoeffs_raw_inactive=allrois_corrcoef_inactive(allrois_corrcoef_inactive~=1); %average correlation coefficient (r, lag0) across all pairs, without autocorrelation
avg_corrcoef_raw_inactive=mean(allcoeffs_raw);

%plot correlation between INACTIVE ROI pairs
figure(17)
inactive_traces_pearson = (allrois_traces_inactive(1:end,index_totake:index_totake+one_minute))';
inactive_detrended = detrend(inactive_traces_pearson,2,round(length(inactive_traces_pearson)*[0.25 0.5 0.75 1]));
[R_inactive,PValue] = corrplot(inactive_detrended);%,'tail','right');
fig_corrcoeff_inactive=figure(18);
colormap(redblue)
imagesc(R_inactive)
caxis([-1 1])
colorbar
ylabel ('roi number')
xlabel ('roi number')
title('corr coeff on INACTIVE roi pairs')

%save all necessary info to calculate (below) average corr coeff across
%active ROIs and across inactive ROIs in each FOV, and then plot together results from all FOVs.
rois_tosave=[savepath '/' mouseID ' _' img_root ' _' img_file_tosave '_active_rois.mat'];
save(rois_tosave, 'one_minute', 'index_totake', 'frate','avg_corrcoef_raw')

rois_tosave_inactive=[savepath_inactive '/' mouseID ' _' img_root ' _' img_file_tosave '_inactive_rois.mat'];
save(rois_tosave_inactive, 'one_minute', 'index_totake', 'frate','avg_corrcoef_raw_inactive')

%% plot mean corr coeffs between active ROIs across FOVs

clear all
close all

root_oxlight="C:\Users\oxlight\active_rois\";
root_ctrl="C:\Users\ctrl\active_rois\";

rois_files_oxlight = dir(fullfile(root_oxlight + "*_active_rois.mat"));
rois_files_ctrl = dir(fullfile(root_ctrl + "*_active_rois.mat"));


colorstouse=colormap(jet(length(rois_files_oxlight)));

allfovs_avg_corrcoeff_raw_oxlight=[];
for ii = 1:length(rois_files_oxlight)
    data = load(rois_files_oxlight(ii).name);
    allfovs_avg_corrcoeff_raw_oxlight=[allfovs_avg_corrcoeff_raw_oxlight data.avg_corrcoef_raw];
end
mean_corrcoeff_raw_oxlight=mean(allfovs_avg_corrcoeff_raw_oxlight);

clearvars -except root_ctrl rois_files_ctrl colorstouse allfovs_avg_corrcoeff_raw_oxlight mean_corrcoeff_raw_oxlight

allfovs_avg_corrcoeff_raw_ctrl=[];
for ii = 1:length(rois_files_ctrl)
    clear data
    data = load(rois_files_ctrl(ii).name);
    
    allfovs_avg_corrcoeff_raw_ctrl=[allfovs_avg_corrcoeff_raw_ctrl data.avg_corrcoef_raw];
end
mean_corrcoeff_raw_ctrl=mean(allfovs_avg_corrcoeff_raw_ctrl);

figure(13)
errorbar(1,mean_corrcoeff_raw_oxlight, std(allfovs_avg_corrcoeff_raw_oxlight)./sqrt(length(allfovs_avg_corrcoeff_raw_oxlight)),'Linewidth', 2.5, 'Color', [0 0 0])
hold on
errorbar(2,mean_corrcoeff_raw_ctrl, std(allfovs_avg_corrcoeff_raw_ctrl)./sqrt(length(allfovs_avg_corrcoeff_raw_ctrl)),'Linewidth', 2.5, 'Color', [0 0 0])
hist_data = [mean_corrcoeff_raw_oxlight mean_corrcoeff_raw_ctrl];
hold on
bar(hist_data, 'k')
hold on
scatter((ones(1,length(allfovs_avg_corrcoeff_raw_oxlight))).*1.25, [allfovs_avg_corrcoeff_raw_oxlight], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [.7 .7 .7])
hold on
scatter((ones(1,length(allfovs_avg_corrcoeff_raw_ctrl))).*2.25, [allfovs_avg_corrcoeff_raw_ctrl], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [.7 .7 .7])
%ylim([-0.2 0.6])
xlim([0.5 2.5])
title('mean corr coef between raw rois during minute of highest activity')
ylabel('mean corr coef')

[h,p] = ttest2(allfovs_avg_corrcoeff_raw_oxlight, allfovs_avg_corrcoeff_raw_ctrl)
[p,h] = ranksum(allfovs_avg_corrcoeff_raw_oxlight, allfovs_avg_corrcoeff_raw_ctrl)

%% compare  mean corr coeffs in oxlight active rois and oxlight inactive rois

oxlight_awake_goodrois=[0.2481    0.3846    0.2868    0.2667    0.2751    0.3230    0.3904    0.1783];
oxlight_awake_darkrois=[0.2162    0.2615    0.2595    0.2161    0.2423    0.2507    0.2364    0.0991];
figure(100)
errorbar(2,mean(oxlight_awake_darkrois), std(oxlight_awake_darkrois)./sqrt(length(oxlight_awake_darkrois)),'Linewidth', 2.5, 'Color', [0 0 0])
hold on
errorbar(1,mean(oxlight_awake_goodrois), std(oxlight_awake_goodrois)./sqrt(length(oxlight_awake_goodrois)),'Linewidth', 2.5, 'Color', [0 0 0])
hist_data = [mean(oxlight_awake_goodrois) mean(oxlight_awake_darkrois)];
hold on
bar(hist_data, 'k')
hold on
scatter((ones(1,length(oxlight_awake_darkrois))).*2.25, [oxlight_awake_darkrois], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [.7 .7 .7])
hold on
scatter((ones(1,length(oxlight_awake_goodrois))).*1.25, [oxlight_awake_goodrois], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [.7 .7 .7])
xlim([0.5 2.5])
title('mean corr coeff in wakefulness on oxlight dark and good rois')

[h,p] = ttest2(oxlight_awake_darkrois, oxlight_awake_goodrois)
[p,h] = ranksum(oxlight_awake_darkrois, oxlight_awake_goodrois)

