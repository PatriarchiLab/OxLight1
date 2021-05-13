%list all the data to analyse. They should be organized as:
% 1) "mouseID": the ID of the mouse
% 2) "image_root": date of the experiment
% 3) "img_file_a": name of the t-series to analyse
% 4) "img_file_tosave": how the t-series will be renamed before saving new data

function [img_root, img_file_a, img_file_tosave, mouseID]=twoP_FOVs_used(oxlight, FOV_toanalyse)

if oxlight == 1
    
    if FOV_toanalyse == 1
        img_root='20201211';
        img_file_a='area01_awakeToAsleep-6629';
        img_file_tosave='area01_asleeptoawake';
        mouseID='073-6220';
    elseif FOV_toanalyse == 2
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 3
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 4
        img_root='';
        img_file_a=''; 
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 5
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 6
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 7
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';        
    elseif FOV_toanalyse == 8
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
        
    end
    
elseif oxlight == 0
    
    if FOV_toanalyse == 1
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 2
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 3
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 4
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 5
        img_root='';
        img_file_a=''; 
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 6
        img_root='';
        img_file_a=''; 
        img_file_tosave='';
        mouseID='';
    elseif FOV_toanalyse == 7
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';        
    elseif FOV_toanalyse == 8
        img_root='';
        img_file_a='';
        img_file_tosave='';
        mouseID='';        
    end
end

end

