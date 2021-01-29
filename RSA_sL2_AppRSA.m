%% RSA Searchlight 
% 2. apply rsa and compute rank correation
% by SPARK 10.Oct.2018

close all
clear all
clc

%% Setup
[ProjSet, fs, Nses, ROIs, Nperm, fname]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.rsatoolbox);
addpath(ProjSet.func);
[subjL, subjN]=callsubj27;
BrRDM='Mtv_DMa24';
UsrOption.BhPath=[ProjSet.Respath, BrRDM, fs];
SvPath=UsrOption.BhPath; %save path
savename=fname.savename; %'sl_RSA_tauA.mat';
if ~exist(SvPath)
    mkdir(SvPath);
end
UsrOption.BrPath=[UsrOption.BhPath, 'Vols', fs];

template24=rdmset(24);
% Pick which model to test
% [E_rdm, Ctx_rdm, Gr_rdm, D_rdm, I_rdm, E_rdm_perm, D_rdm_perm, I_rdm_perm]=EucRDM_factorize_ss
[modelRDM.org, ~, ~, ~, ~, modelRDM.perm, ~, ~]=EucRDM_factorize_ss(template24, Nperm); modeltype='Euc'; % Euclidean
% [~, ~, ~, modelRDM.org, ~, ~, modelRDM.perm, ~]=EucRDM_factorize_ss(template24, Nperm); modeltype='D'; % 1-D rank distance in relevant dimension
% [~, ~, ~, ~, modelRDM.org, ~, ~, modelRDM.perm]=EucRDM_factorize_ss(template24, Nperm); modeltype='I'; % 1-D rank distance in irrelevant dimension
% call EucRDM_factorize_pc instead of EucRDM_factorize_ss for partial-correalation out 
% e.g.[modelRDM.org, ~, ~, ~, ~, modelRDM.perm, ~, ~]=EucRDM_factorize_pc(template24, Nperm); modeltype='Euc'; % Euclidean
modelRDMs=[{modelRDM.org}; modelRDM.perm]';

%% Turn on the parallel computing
pflag = gcp('nocreate');
if isempty(pflag)
    poolsize=0;
else
    poolsize=pflag.NumWorkers;
end

%% Main 
BrSUB = dir(fullfile(UsrOption.BrPath,'*_sL_RDM.mat')); % results from RSA_sL1_RDMcrs.m    
    for sbi=1:size(BrSUB,1)
        clear curSub
        curSub=BrSUB(sbi).name;
        [kendallT] = slctSubModel(curSub, modelRDMs, UsrOption);
        T0.(modeltype)(sbi,:)=kendallT(1,:); % Testing model
        B.(modeltype)(sbi,:)=mean(kendallT(2:Nperm+1,:)); % Suffled baseline
        tau.(modeltype)(sbi,:)=T0.(modeltype)(sbi,:)-B.(modeltype)(sbi,:);
    end
tau_idx{1}='Kendalls tauA';
save([SvPath, savename], 'tau', 'tau_idx'); 
cd(SvPath);
delete(gcp('nocreate'));

%% custom functions %
function [CorrT] = slctSubModel(curSub, modelRDMs_cell_orth, UsrOption)  
    parfor mdi=1:size(modelRDMs_cell_orth,2)     % size of NPerm+1                                   
        curRDM=modelRDMs_cell_orth{mdi};    %current RDM
        [taua] = appRSA(curSub, curRDM, UsrOption);
        CorrT(mdi,:) = taua';
    end
end

function [taua] = appRSA(sub, rdm, usroption)
    curROI=load([usroption.BrPath, sub]); %RDMs
    RDMs=curROI.RDMs;
    parfor vi=1:size(RDMs,2) % number of searchlights                                                                                      
        taua(vi)=rsa.stat.rankCorr_Kendall_taua(selectriu(RDMs(vi).RDMeuc), selectriu(rdm)); 
    end
end


