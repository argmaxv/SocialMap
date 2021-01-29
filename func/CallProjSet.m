
function [ProjSet, fs, Nses, ROIs, Nperm, fname]=CallProjSet
    
    fs=filesep;
    Nses=2;
    Nperm=1000;% number of permutation
    ROIs={'HCr', 'HCl', 'ECr', 'ECl', 'mOFCr', 'mOFCl', 'AMr', 'AMl', 'M1r', 'M1l'};
    
    ProjSet.basepath='/mnt/datashare/Pr3.SS/Imaging';
    ProjSet.spmdir='/usr/local/MATLAB/spm12';
    ProjSet.rsatoolbox='/mnt/datashare/Pr3.SS/Programs/rsatoolbox-develop';
    ProjSet.func=[ProjSet.basepath, fs, 'Batch', fs, 'func']; % where custom functions are
    ProjSet.ONSETpath=[ProjSet.basepath, fs, 'Onset', fs]; %behavioral data
    ProjSet.DATApath=[ProjSet.basepath, fs, 'Data', fs];    %raw data
    ProjSet.ANApath=[ProjSet.basepath, fs, 'Analysis', fs]; 
    ProjSet.ANA1path=[ProjSet.ANApath, fs, 'Analysis Lv1', fs]; %lv1 analysis
    ProjSet.ANA2path=[ProjSet.ANApath, fs, 'Analysis Lv2', fs]; %lv2 analysis
    ProjSet.ROIpath=[ProjSet.ANApath, 'Analysis ROI', fs];  %ROI path
    ProjSet.Maskpath=[ProjSet.ROIpath, 'Mask', fs]; 
    ProjSet.Respath=[ProjSet.ANApath, 'Analysis RSA', fs]; % save folder for RSA searchlight
    
    fname.savename='sl_RSA_tauA.mat';
    fname.savefolder='TauMap_Euc';
    fname.maskname='Mask.img';
    fname.searchlightname='searchlight_100.mat';
    
end




