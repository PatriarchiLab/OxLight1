% use this script to:
% 1)upload, for each FOV, data saved using the script binarize_deconvolve
% 2)plot standard deviation changes during asleep and awakening periods, across FOVs
%3) plot F0 for oxlight and ctrl FOVs.
%***these plots are persent in extended figure 8, panels E and F of Duffet et al.

%% upload data
clear all
close all

root_oxlight="C:\Users\oxlight\deconvolution\";
root_ctrl="C:\Users\ctrl\deconvolution\";

binarization_files_oxlight = dir(fullfile(root_oxlight + "*_deconvolved.mat"));
binarization_files_ctrl = dir(fullfile(root_ctrl + "*_deconvolved.mat"));

%% run for oxlight and control mice separately

avg_std_asleep=[];
avg_std_awake=[];
all_F0_oxlight=[];
for i = 1:length(binarization_files_oxlight)
    clear data
    data = load(binarization_files_oxlight(i).name);
    
    avg_std_asleep=[avg_std_asleep  mean(data.std_toplot(1:floor(length(data.std_toplot)/2)))];
    avg_std_awake=[avg_std_awake  mean(data.std_toplot(floor(length(data.std_toplot)/2):end))];
    
    figure(1); hold on
    scatter(1,data.F0, 'g', 'filled')
    all_F0_oxlight=[all_F0_oxlight data.F0];
end

clearvars -except all_F0_oxlight avg_std_asleep avg_std_awake root_ctrl binarization_files_ctrl

avg_std_asleep_ctrl=[];
avg_std_awake_ctrl=[];
all_F0_ctrl=[];
for i = 1:length(binarization_files_ctrl)
    clear data
    data = load(binarization_files_ctrl(i).name);
    binarization_files_ctrl(i).name
    
    avg_std_asleep_ctrl=[avg_std_asleep_ctrl  mean(data.std_toplot(1:floor(length(data.std_toplot)/2)))];
    avg_std_awake_ctrl=[avg_std_awake_ctrl  mean(data.std_toplot(floor(length(data.std_toplot)/2):end))];
    
    figure(1); hold on
    scatter(1.3,data.F0, 'k', 'filled')
    xlim([0.9 1.4])
    xticks([1 1.3])
    set(gca,'TickDir', 'out', 'FontSize', 14, 'Linewidth', 2.5, 'xticklabels', {'oxlight' 'ctrl'})
    ylabel('fluorescence (a.u.)')
    all_F0_ctrl=[all_F0_ctrl data.F0];
    
end

%stats on F0 data
[p_F0,h_F0] = ttest(all_F0_oxlight, all_F0_ctrl)

%plot std changes and compare them
figure (10)
plot([1 1.3],[mean(avg_std_asleep) mean(avg_std_awake)], '-', 'Color', [0 1 0], 'LineWidth', 2.5)
hold on
plot([1 1.3],[mean(avg_std_asleep_ctrl) mean(avg_std_awake_ctrl)], '-', 'Color', [0 0 0], 'LineWidth', 2.5)
box off
hold on
for ww= 1:length(avg_std_asleep)
    plot([1 1.3],[avg_std_asleep(ww) avg_std_awake(ww)], '--', 'Color', [0 1 0], 'LineWidth', 1)
end
hold on
for hh= 1:length(avg_std_asleep_ctrl)
    plot([1 1.3],[avg_std_asleep_ctrl(hh) avg_std_awake_ctrl(hh)], '--', 'Color', [0 0 0], 'LineWidth', 1)
end
ylim([0 200])
xlim([0.95 1.35])
ylabel ('mean sd')
set(gca,'TickDir', 'out', 'FontSize', 14, 'Linewidth', 2.5, 'xticklabels', {'asleep', 'awake'})
xticks([1 1.3])
% stats
[p_oxlight,h_oxlight] = ttest(avg_std_asleep, avg_std_awake)
[p_ctrl,h_ctrl] = ttest(avg_std_asleep_ctrl, avg_std_awake_ctrl)
[pp_oxlight,hh_oxlight] = ranksum(avg_std_asleep, avg_std_awake)
[pp_ctrl,hh_ctrl] = ranksum(avg_std_asleep_ctrl, avg_std_awake_ctrl)

