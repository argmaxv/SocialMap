%% RSA Searchlight 
% 3. rankCorrelation to the map and smoothing.
% by SPARK 10.Oct.2018
%
% read KendallT (T.(RDM) has sub x voxels matrix)
% write to nii images (MaskMap size)
% smooth individual tau map

close all
clear all
clc

%% Setup
[ProjSet, fs, Nses, ROIs, Nperm, fname]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.rsatoolbox);
addpath(ProjSet.func);
BrDataPath=ProjSet.Respath;
ModelNames={'Mtv_DMa24'};
dataname=fname.savename; % sl_RSA_tauA.mat;
svfolder=fname.savefolder;     % TauMap_Euc;

% Mask info to make 1D array to 3D map
Mask = fullfile(ProjSet.Maskpath, fname.maskname);
load(fullfile(ProjSet.Maskpath, fname.searchlightname)); %voxel
struct=spm_vol(deblank(Mask));
MaskMap =spm_read_vols(struct);
vsize=struct.dim;

nVox=length(voxel);
FWHMmm=8;
prefix='s';

%% turn on the parallel computing
pflag = gcp('nocreate');
if isempty(pflag)
    poolsize=0;
else
    poolsize=pflag.NumWorkers;
end
    
%% Main
for mm=1:numel(ModelNames)
    ModelName=ModelNames{mm};
    datapath=[BrDataPath, ModelName, fs];             % data path
    savepath=[BrDataPath, ModelName, fs, 'sL', fs]; % save path
    load([datapath, dataname, '.mat']); %tau.Euc
    T=tau;
    RDMnames=fieldnames(T); %Euc
    nRDM=size(RDMnames,1);
    clear tau

    for iRDM=1:nRDM %nRDM=1, Euc
        sSimilarity=[];
        disp(RDMnames{iRDM});
        mappath=fullfile([savepath, svfolder], RDMnames{iRDM}); % save path
        if ~exist(mappath)
            mkdir(mappath);
        end
        parfor sub=1:size(T.(RDMnames{iRDM}),1) %%%% ### par
            thearray = T.(RDMnames{iRDM})(sub,:);
            filename = [ModelName, '_', num2str(sub),'_', RDMnames{iRDM}];
            [~, fname] = array2nii(thearray, voxel, MaskMap, struct, mappath, filename); % save unsmoothed nii
            smoothingnii(FWHMmm, prefix, fname, 1); % save smoothed nii with the prefix
        end
    end
    disp('All maps are generated and smoothed.');
end
