% use this script for:
% 1)binarizing each time series according to the thereshold explained in the methods section of Duffet et al.
% 2)deconvolving activity on the binarized time-series
% **the following 3 sections should be run on each individual FOV.

%% set up all parameters and variables before running the script
close all
clear all

%choose the FOV to analyse
oxlight=1; %set to 0 if it is an oxlight control, set to 1 if it is an oxlight mouse
FOV_toanalyse= 7;
[img_root, img_file_a, img_file_tosave, mouseID]=twoP_FOVs_used(oxlight,FOV_toanalyse);

%call the suite2p output file in order to get the frame rate of the time series
load(['C:\Users\suite2p\plane0\Fall.mat'])
frate=ops.fs;
%get the whole time series (registered tiff file) loaded and ready for binarization
root_FOV="C:\Users\suite2p\plane0\reg_tif\";
tiff_folder = dir(fullfile(root_FOV + "*.tif"));

%set the correct parameters
asleeptoawake=1; %1 if mouse starts asleep and ends awake, 0 if mouse starts awake and ends asleep
pixel_size = 1.274;
object_size_um = 50; %size of the localized activity in um
object_px = floor(object_size_um/pixel_size); %size of the localized activity in pixels

if oxlight
    savepath=(['C:\Users\oxlight\deconvolution']);
else
    savepath=(['C:\Users\ctrl\deconvolution']);
end


%% run binarization

full_matrix=[]; %initialize 3d matrix that will contain the whole acquisition
for f = 1:length(tiff_folder)
    hInfo = imfinfo(root_FOV + "/" + tiff_folder(f).name);
    
    %// Set parameters.
    ImageHeight = hInfo(f).Height;
    ImageWidth = hInfo(f).Width;
    SliceNumber = numel(hInfo);
    
    %// Open Tiff object
    Stack_TiffObject = Tiff(root_FOV + "/" + tiff_folder(f).name,'r');
    
    %// Initialize array containing your images.
    ImageMatrix = zeros(ImageHeight,ImageWidth,SliceNumber,'uint32');
    for k = 1:SliceNumber
        
        %// Loop through each image
        Stack_TiffObject.setDirectory(k)
        
        %// Put it in the array
        ImageMatrix(:,:,k) = Stack_TiffObject.read();
    end
    
    full_matrix=cat(3,full_matrix, ImageMatrix);
    
    %// Close the Tiff object
    Stack_TiffObject.close
end

%calculate baseline fluorescence on the time series.
avg_trace_fullmatrix = squeeze(mean(mean(full_matrix, 1), 2))';
figure (1)
subplot(2,1,1)
plot(avg_trace_fullmatrix)
[transition, dFoF,F0, F0_std]=twoP_baseline_bistate(avg_trace_fullmatrix,frate,asleeptoawake);
subplot(2,1,2)
plot(dFoF)

%calculate standard deviation across frames during sleeping period (first minute of t-series) and more active periods (eg. most active minute or last minute of imaging)
secs_long=length(dFoF)/frate; %how long is the acquisition in seconds
one_minute=floor((length(dFoF)*60)/secs_long); %how many frames correspond to one minute of the series.
sleeping_period=full_matrix(:,:,1:one_minute);
allframes_sleeping_std=[];
for hh=1:size(sleeping_period,3)
    sleeping_std_framebyframe=std(double(sleeping_period(:,:,hh)),0, 'all');
    allframes_sleeping_std=[allframes_sleeping_std sleeping_std_framebyframe];
end
non_sleeping_period=full_matrix(:,:,end-(one_minute):end);
allframes_std=[];
for kk=1:size(non_sleeping_period,3)
    std_framebyframe=std(double(non_sleeping_period(:,:,kk)),0, 'all');
    allframes_std=[allframes_std std_framebyframe];
end
std_toplot=[allframes_sleeping_std allframes_std];
[h_value_std,p_value_std] = ttest2(allframes_sleeping_std, allframes_std, 'tail', 'left');

%binarize part of the time series
secs_long=length(dFoF)/frate; %how long is the acquisition in seconds
one_minute=floor((length(dFoF)*60)/secs_long); %how many frames correspond to one minute of the series.
period_to_binarize=one_minute*2:one_minute*6; %choose the range of frames to binarize: currently binarization starts one minute after tuning isoflu off, and continues for 4 minutes
matrix_to_binarize=full_matrix(:,:,period_to_binarize); %portion of time series to binarize
avg_original_FOV=mean(matrix_to_binarize,3); %average over frames (ie over time)
max_original_FOV=max(matrix_to_binarize,[],3);
std_original_FOV=std(double(matrix_to_binarize),0,3);
threshold_toplot=mean2(std_original_FOV); %threshold for binarization
figure(2)
imagesc(avg_original_FOV)
colormap(gray)
colorbar
%now binarize each frame
bw_framebyframe_FOV_pxbypx = matrix_to_binarize > avg_original_FOV + threshold_toplot; % every pixel is individually binarized based on its own mean
fname1 = [mouseID '-' img_root '-' img_file_tosave '-one_std_deconv_l1.tif']; %name of the tiff file that we're going to write
fname_blurred = [mouseID '-' img_root '-' img_file_tosave '-one_std_blurred.tif'];

%% run deconvolution

PSF_1 = fspecial('gaussian',[object_px object_px],object_px/2);
PSF_smooth = fspecial('gaussian',[size(matrix_to_binarize,1) size(matrix_to_binarize,2)],2);
figure(3),imshow(PSF_1,[])

%write tiff file from the binarized matrix
i_img_1 = [];
cmp = gray;
iwrite = true;
parfor nn = 1:size(matrix_to_binarize, 3)
    %deconvolve each frame
    blurred(:,:,nn) = imfilter(double(bw_framebyframe_FOV_pxbypx(:,:,nn)),PSF_smooth,'symmetric','conv'); %apply gausian filter for smoothing
    I = edgetaper(double(blurred(:,:,nn)),PSF_1);
    luc1(:,:,nn) = deconvlucy(I,PSF_1,5); %only run 5 iterations, because non-linear deblurring deteriorates after too many iterations
end

max_luc1 = max(luc1,[], 'all');
min_luc1 = min(luc1,[], 'all');
max_blurred = max(blurred,[], 'all');
min_blurred = min(blurred,[], 'all');
avg_blurred=mean(blurred, 3);
avg_luc1=mean(luc1, 3);
[mxv,idx] = max(luc1(:));
[r,c,p] = ind2sub(size(luc1),idx);
for hh = 1:size(luc1, 3)
    norm_luc1(:,:,hh) = (luc1(:,:,hh)-min_luc1)./(max_luc1-min_luc1);
    norm_blurred(:,:,hh) = (blurred(:,:,hh)-min_blurred)./(max_blurred-min_blurred);
end
std_luc1=std(luc1(:,:,p),[], 'all');
figure(4); imagesc(blurred(:,:,p))
title('frame with highest activity - blurred')
colorbar
figure(5); imagesc(luc1(:,:,p))
title('frame with highest activity - deconvolved')
colorbar
figure(6); imagesc(avg_blurred)
title('averaged blurred')
colorbar

%visualize the average projection over frames (ie over time) onf the deconvolved FOV. Also, check distribution of pixel values in the projection, and in all individual frames.
img_avg_luc1=figure(7); imagesc(avg_luc1)
title('averaged deconvolved')
colorbar
figure(8) %hist on just the average projection
edges = [min(avg_luc1,[],'all'):0.01:max(avg_luc1,[],'all')];
h_avg = histogram(avg_luc1, edges);
clear edges h
img_allpixelsvalues=figure(9); %hist all pixels from all frames (together)
edges = [min_luc1:0.05:max_luc1];
h_allframes = histogram(luc1, edges);

%count (in percentage) how many pixels on the deconvolved average
%projection have a value higher than either mean+sd or mean+2sd
matrix_aboveThresh_1sd = avg_luc1 > mean2(avg_luc1)+1*std2(avg_luc1);
perc_aboveThresh_1sd = sum(matrix_aboveThresh_1sd, 'all')*100/(size(avg_luc1,1)*size(avg_luc1,2));
matrix_aboveThresh_2sd = avg_luc1 > mean2(avg_luc1)+2*std2(avg_luc1);
perc_aboveThresh_2sd = sum(matrix_aboveThresh_2sd, 'all')*100/(size(avg_luc1,1)*size(avg_luc1,2));
img_objects_abovethresh=figure(10);
subplot(1,2,1)
imagesc(matrix_aboveThresh_1sd)
title('mean + 1sd')
subplot(1,2,2)
imagesc(matrix_aboveThresh_2sd)
title('mean + 2sd')

%save everything
deconvolved_tosave=[savepath '/' mouseID ' _' img_root ' _' img_file_tosave '_deconvolved.mat'];
save(deconvolved_tosave,'frate', 'F0', 'one_minute', 'period_to_binarize', 'luc1', 'blurred', 'pixel_size', 'object_size_um', 'object_px', 'avg_luc1', 'h_avg', 'h_allframes', 'matrix_aboveThresh_1sd', 'perc_aboveThresh_1sd', 'matrix_aboveThresh_2sd', 'perc_aboveThresh_2sd', 'std_toplot')

%% for each individual FOV, you can run this section to create a video of the deconvolved time-series

for gg = 1:size(luc1, 3)
    luc1_rgb_video(:,:,:,gg) = ind2rgb(round(norm_luc1(:,:,gg)*256), cmp);
    blurred_rgb_video(:,:,:,gg) = ind2rgb(round(norm_blurred(:,:,gg)*256), cmp);
    
    % Make an RGB image:
    i_img_1 = ind2rgb(round(norm_luc1(:,:,gg)*256), cmp);
    i_img_blurred = ind2rgb(round(norm_blurred(:,:,gg)*256), cmp);
    
    % Generate your tiff stack:
    if iwrite
        if gg == 1
            % First slice:
            imwrite(i_img_1,fname1);
            imwrite(i_img_blurred,fname_blurred);
        else
            % Subsequent slices:
            imwrite(i_img_1,fname1,'WriteMode','append');
            imwrite(i_img_blurred,fname_blurred,'WriteMode','append');
        end
    end
    
    disp(gg)
    
end
implay(blurred_rgb_video,80);
implay(luc1_rgb_video,80);

