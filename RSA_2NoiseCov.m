%% RSA ROI
% 2. NoiseNormalization and generate RDM per ROI
% by SPARK 1.Oct.2018

clear all;
close all;
clc;

%% Setting
[ProjSet, fs, Nses, ROIs]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.rsatoolbox);
addpath(ProjSet.func);
raw_path=ProjSet.DATApath; %Data path
Anal1path=ProjSet.ANA1path; % Lv1 analysis spm.mat
[subj, nsubj]=callsubj27;
    
models={'Mtv_DMa24'};
colors={[1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [0,1,1], [1,1,1]};
nEtc=7; % to count the number of betas 7=6: motion regressors + 1: session mean;

%% Main
for m=1:numel(models)
    modelname=models{m};
        
    for roiIdx=1:numel(ROIs) % Loop ROI
        
        ROIselec=roiIdx;
        ROIpath=ROIs{roiIdx};
        svpath=[raw_path, ROIpath, fs, 'noiseOut', fs, [modelname, '_', datestr(date)], fs]; %save path
        if ~exist(svpath, 'dir')
            mkdir(svpath);
        end    
        disp(['*** the current ROI is ', ROIpath, ' ***']);
        clear RDMs

        for s=1:nsubj
            clear rawdata u_hat nBeta zUhat Uhat
            rawdata=load([raw_path, ROIpath, fs, subj{s}, '_raw.mat']); %from the outcomes of RSA_1ExtctRaw
            Y=rawdata.RAW.Traw;
            cor.xyz=rawdata.RAW.xyz;
            cor.mni=rawdata.RAW.mni;
            clear SPM spmpath design
            spmpath=[Anal1path, modelname, fs, subj{s}, fs, 'SPM.mat'];
            design=load(spmpath);
            SPM=design.SPM;
            [u_hat,resMS,Sw_hat,beta_hat,shrinkage,trRR]=rsa.spm.noiseNormalizeBeta(Y,SPM); % RSA toolbox
            u_hat=real(u_hat);
            nBeta=size(u_hat,1); % compute mean u_hat across session
            nBeta=(nBeta - (Nses*nEtc))/Nses;
            
            for nR=1:Nses
                stp = 1 + (nR-1)*(nBeta+nEtc-1);
                enp = stp + nBeta -1;
                Uhat{nR}=u_hat(stp:enp,:);
                for pnt=stp:enp
                    zUhat{nR}(pnt-(stp-1),:) = (u_hat(pnt,:)- nanmean(u_hat(pnt,:)) ) / nanstd(u_hat(pnt,:));
                end

% Within Session RDM 
                rdm_zses{nR}=squareform(pdist(zUhat{nR},'Euclidean'));
                RDMs(1,s,nR).RDM=rdm_zses{nR};
                RDMs(1,s,nR).name=[models{m},'_', ROIpath, ' | ' subj{s}, ' | Session: ' num2str(nR)];
                RDMs(1,s,nR).color=colors{m};
            end %for nR
            meanRDMwth=rsa.rdm.averageRDMs_subjectSession(RDMs(1,s), 'session');
            RDMsanity(s).EDIwth = extractEDIs(meanRDMwth.RDM);
            RDMsanity(s).RDMwth=meanRDMwth.RDM; 
            RDMwth(1,s).RDM=meanRDMwth.RDM;
            RDMwth(1,s).name=[models{m},'_', ROIpath, ' | ' subj{s}, ' | Session: 1'];
            RDMwth(1,s).color=colors{m};
            
% Cross Session RDM
            clear crsSesRDM_ses crsSesRDM;
            rns=[1:Nses];
            crsSesRDM_sum=zeros(size(zUhat{1},1));
            for r0=1:Nses
                rns_cur=rns(rns~=r0);
                for r0_cur=1:numel(rns_cur)
                    crsSesRDM_ses{r0}=pdist2(zUhat{rns(r0)},zUhat{rns_cur(r0_cur)},'Euclidean');
                end
                crsSesRDM_sum=crsSesRDM_ses{r0}+crsSesRDM_sum;
            end
            crsSesRDM=crsSesRDM_sum/Nses;
            RDMsanity(s).EDIcrs = extractEDIs(crsSesRDM);          
            RDMsanity(s).RDMcrs=crsSesRDM;
            RDMcrs(1,s).RDM=crsSesRDM;
            RDMcrs(1,s).name=[models{m},'_', ROIpath, ' | ' subj{s}, ' | Session: 1'];
            RDMcrs(1,s).color=colors{m};                

            disp([subj{s}, ' multivariate noise normalization done']);
            clear Uhat u_hat resMS Sw_hat beta_hat shrinkage trRR rawdata Y spmpath design SPM cor;

        end %Subject
        save([svpath, 'BrainRDM.mat'], 'RDMs', 'RDMwth', 'RDMcrs', 'RDMsanity');
    end %for ROIs
end % for Models
