% Use this script to compute a distribution of pixels values across all deconvolved FOVs.
% Identify as belonging to "active ROIs" the regions where pixel values are between the 80th and the 100th of the distribution
% Identify as belonging to "inactive ROIs" the regions where pixel values are between the 6th and the 69th of the distribution
% all ROIs will be defined manually on the median projection of the FOV
% visualised in suite2p. In Duffet et al, first 10 ROIs identified are the inactive ones, last 10 identified are the active ones.
% **the files used in this script are those saved in the binarize_deconvolve script

%%
clear all
close all

root_oxlight="C:\Users\oxlight\deconvolution\";
root_ctrl="C:\Users\ctrl\deconvolution\";

deconv_files_oxlight = dir(fullfile(root_oxlight + "*_deconvolved.mat"));
deconv_files_ctrl = dir(fullfile(root_ctrl + "*_deconvolved.mat"));

colorstouse=colormap(jet(length(deconv_files_oxlight)));

%pixel distributions for all oxlight deconvolved FOVs
alloxlight_distrib_inprojection=[];
for ii = 1:length(deconv_files_oxlight)
    clear data
    data = load(deconv_files_oxlight(ii).name);
    
    %histogram all pixels across all frames
    img_allpxvalues_oxlight=figure(1); hold on
    plot(data.h_allframes.BinEdges(2:end), data.h_allframes.BinCounts, 'Linewidth', 0.5, 'Color', colorstouse(ii, :))
    
    alloxlight_distrib_inprojection=[alloxlight_distrib_inprojection reshape(data.avg_luc1,1, size(data.avg_luc1,1)*size(data.avg_luc1,2))];
end
figure(1)
title('oxlight - deconv activity across frames')
xlabel('deconvolved activity')
ylabel('counts')
xlim([0 2])
ylim([0 20000000])

clearvars -except root_ctrl deconv_files_ctrl colorstouse alloxlight_distrib_inprojection

%pixel distributions for all ctrl deconvolved FOVs
allctrl_distrib_inprojection=[];
for ii = 1:length(deconv_files_ctrl)
    clear data
    data = load(deconv_files_ctrl(ii).name);
    
    img_allpxvalues_ctrl=figure(10); hold on
    plot(data.h_allframes.BinEdges(2:end), data.h_allframes.BinCounts, 'Linewidth', 0.5, 'Color', colorstouse(ii, :))
    
    allctrl_distrib_inprojection=[allctrl_distrib_inprojection reshape(data.avg_luc1,1, size(data.avg_luc1,1)*size(data.avg_luc1,2))];
end
figure(10)
title('ctrl - deconv activity across frames')
xlabel('deconvolved activity')
ylabel('counts')
xlim([0 2])
ylim([0 20000000])

%%%% set thresholds for defining ROIs
distrib_allFOVs=[allctrl_distrib_inprojection alloxlight_distrib_inprojection];
highpass_active_rois=prctile(distrib_allFOVs,80); %active rois are drawn where the average deconvolved projection has this value or above
bandpass_nonactiverois=prctile(distrib_allFOVs,[6 79]); %non active rois are taken between these two values
%%%%
