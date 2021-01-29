%% RSA Searchlight 
% 0. Set n serchlight including 100 voxels
% by SPARK 10.Oct.2018

%% Setting
% The path where the mask (Mask.img) is
[ProjSet, fs, Nses, ROIs, Nperm, fname]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.rsatoolbox);
maskpath=ProjSet.Maskpath;
svpath=maskpath; %save path 

%% Main
Vmask=spm_vol(fullfile(maskpath, fname.maskname));
Vmask.data=spm_read_vols(Vmask);
L = rsa.defineSearchlight({Vmask},Vmask,'shere',[15 100]);
save(fullfile(svpath, fname.searchlightname) ,'-struct', 'L' );
