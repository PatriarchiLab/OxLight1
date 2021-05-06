
% Import the configuration options
conf = DFFanalysisConfig();

% Add necessary paths
% Warning: call addpath only once because it is slow!
paths = strcat(conf.DFFanalysisPath, ';', ...
               conf.DFFanalysisPath,'/core;', ...
               conf.DFFanalysisPath,'/core/bfmatlab;');
addpath(paths);

% Run FRETanalysis
fDFFanalysis(conf);