function run_kilosort_onfile(rootZ,n_channels)
%% Based on code by Marius Pachitariu (marius10p, github), Kilosort2 code
%% Requires Kilosort2 and its dependencies
%
%Joseph M. Stujenske, 2022
%
if nargin<2 || isempty(n_channels)
    n_channels=28;
end
addpath(genpath('C:\Users\Admin\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath('C:\Users\Admin\Documents\MATLAB\npy-matlab-master') % for converting to Phy
% rootZ = 'C:\Users\Admin\Documents\MATLAB\Testtest'; % the raw data binary file is in this folder
rootH = 'C:\Users\Admin\Documents\'; % path to temporary binary file (same size as data, should be on fast SSD)
pathToYourConfigFile = 'C:\Users\Admin\Documents\MATLAB\Kilosort2-master\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'Ch28.mat';
n_channels=28;
%automatically detect noisy channels
M=memmapfile([rootZ,'\data_binary.bin']);
M=memmapfile([rootZ,'\data_binary.bin'],'Format',{'int16',[n_channels length(M.Data)/n_channels/2],'data'});
A=M.Data.data(:,1:min(100000,size(M.Data.data,2)));
ranges=max(A,[],2)-min(A,[],2);
noisy=ranges>min(double(ranges))*4;
load([pathToYourConfigFile,'\Ch28.mat'],'chanMap','chanMap0ind','connected','kcoords','name','xcoords','ycoords');
connected=true(28,1);
connected(noisy)=false;
save([pathToYourConfigFile,'\Ch28.mat'],'chanMap','chanMap0ind','connected','kcoords','name','xcoords','ycoords');
%
ops.trange = [0 Inf]; % time range to sort
ops.NchanTOT    = n_channels; % total number of channels in your recording

run(fullfile(pathToYourConfigFile, 'configFile32.m'))
ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);

%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ);

%% if you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

% Ensure all GPU arrays are transferred to CPU side before saving to .mat
rez_fields = fieldnames(rez);
for i = 1:numel(rez_fields)
    field_name = rez_fields{i};
    if(isa(rez.(field_name), 'gpuArray'))
        rez.(field_name) = gather(rez.(field_name));
    end
end

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(rootZ, 'rez2.mat');
save(fname, 'rez', '-v7.3');
