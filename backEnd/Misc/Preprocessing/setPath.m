%% Test to see if current folder is correct

currentFolder = pwd;
if ~strcmp(currentFolder(end-25:end), 'DiffusionLearning-INTEGRAL')
    error(['You are in the wrong directory.' newline 'Change your current folder to:' newline '"DiffusionLearning-INTEGRAL".'])
end

%% Restore default path
% this gets rid of irrelevant files that may interfere with our code

restoredefaultpath

%% Add folders in directory to path

addpath(genpath(currentFolder))

%% Get rid of old results from path

rmpath(genpath(strcat(currentFolder, '/results/oldResults')))