%% RSA ROI
% 3. Compute Kendall's tauA from the patterns of each ROI
% by SPARK 1.Oct.2018

clear all
close all
clc

%% Setting
[ProjSet, fs, Nses, ROIs, Nperm]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.rsatoolbox);
addpath(ProjSet.func);
datapath = ProjSet.DATApath; %Data path 
resRSApath=ProjSet.ROIpath; %Save path 
[subj_selec, nsubj_selec]=callsubj27;
nROIs=numel(ROIs);
fltrpathlist = {'Mtv_DMa24'};
svoption=0; %1=save the results
rdmtype={'Org','PC'};
OrgPc=1; % to indicate the type of analyis in the results file name

%% Main
for rdm_model=1:numel(fltrpathlist)
    fltrpath = fltrpathlist{rdm_model};
    
% Define Behavior RDM      
        template24=rdmset(24);
        [E_rdm, Ctx_rdm, Gr_rdm, D_rdm, I_rdm, E_rdm_perm, D_rdm_perm, I_rdm_perm]=EucRDM_factorize_ss(template24, Nperm); %Generate RDM + permuted RDM for D,E, and I
        %[E_rdm, Ctx_rdm, Gr_rdm, D_rdm, I_rdm, E_rdm_perm, D_rdm_perm, I_rdm_perm]=EucRDM_factorize_pc(template24, Nperm); %Generate RDM (partial corr out)
        bhv_rdm_idx={'D_rdm','Ctx_rdm','Gr_rdm','E_rdm','I_rdm'};
        for rn=1:numel(bhv_rdm_idx)
            if ~exist([bhv_rdm_idx{rn}, '_perm'], 'var')
                eval([bhv_rdm_idx{rn}, '_perm=RDMshuffles(', bhv_rdm_idx{rn} ',Nperm);']); % generate permuted RDM for Context and Group
            end
            eval(['bhv_rdm{1,rn}=', bhv_rdm_idx{rn}, ';']);
            eval(['bhv_rdm_perm{:,rn}=', bhv_rdm_idx{rn}, '_perm;']);
        end
        bhv_rdm_numreg = size(bhv_rdm,2); % number of regressors

% Initialize variables for collecting stats across ROIs 
        collec.p=[];
        collec.tau_mean=[];
        collec.tau_se=[];

% Define Brain RDM
        for oi=1:numel(ROIs)
                theROI=ROIs{oi};
                rdmpath=[datapath, theROI, fs, 'noiseOut', fs, fltrpath, fs];
                disp([theROI, '     ', fltrpath]);
                load([rdmpath fs 'BrainRDM.mat']); 
                brainRDMs=RDMcrs;

                for bx=1:size(brainRDMs,2) %Loop subjects
                        for mx=1:bhv_rdm_numreg % n model RDMs
                            tau0(bx,mx)=rsa.stat.rankCorr_Kendall_taua(selectriu(brainRDMs(bx).RDM), selectriu(bhv_rdm{1, mx})); %Initial tauA
                            for pm=1:Nperm
                                tau_perm(bx,mx,pm)=rsa.stat.rankCorr_Kendall_taua(selectriu(brainRDMs(bx).RDM), selectriu(bhv_rdm_perm{1, mx}{pm,1})); %Baseline
                            end
                            tau(bx,mx)=tau0(bx,mx)-mean(tau_perm(bx,mx,:)); %Kendall's tauA
                        end
                end

                for mx=1:bhv_rdm_numreg % Group level analysis
                    tau_p(mx)=rsa.stat.signrank_onesided(tau(:,mx)); 
                    tau_mean(mx) = nanmean(tau(:,mx));
                    tau_se(mx) = SEM(tau(:,mx));
                end
% Write results to 'stats'
                stats.(ROIs{oi}).(rdmtype{OrgPc}).tau=tau;
                stats.(ROIs{oi}).(rdmtype{OrgPc}).tau_mean=tau_mean;
                stats.(ROIs{oi}).(rdmtype{OrgPc}).tau_se=tau_se;
                stats.(ROIs{oi}).(rdmtype{OrgPc}).p=tau_p;
                stats.(ROIs{oi}).(rdmtype{OrgPc}).Idx=bhv_rdm_idx;
                collec.p=[collec.p; tau_p];
                collec.tau_mean=[collec.tau_mean; tau_mean];
                collec.tau_se=[collec.tau_se; tau_se];
        end % for ROIs
        
        [~,pFWE]=bonferroni_holm((collec.p(:,1:end-1))); %multiple comparison correction (Models (Cntx,D,E, and Gr) x nROIs)
        stats.All.(rdmtype{OrgPc}).pFWE=pFWE;
        stats.All.(rdmtype{OrgPc}).p=collec.p;
        stats.All.(rdmtype{OrgPc}).tau_mean=collec.tau_mean;
        stats.All.(rdmtype{OrgPc}).tau_se=collec.tau_se;
        
% Save to the file
        if svoption 
            svPath=[resRSApath, fltrpath]; % Save path
            if ~exist(svPath, 'dir')
                mkdir(svPath);
            end
            save(fullfile(svPath, ['RDM_stats', rdmtype{OrgPc}, '.mat']), 'stats'); 
        end
       
end %for type of RSA 
