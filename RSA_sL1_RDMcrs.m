%% RSA Searchlight 
% 1. extract Vols from spmT
% by SPARK 10.Oct.2018
% Create RDM matrix per searchlight
% Dissimilarity is measured by Euclidean distance between activity acquired from differnt sessions 

close all
clear all
clc

%% Setting
[ProjSet, fs, Nses]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.rsatoolbox);
addpath(ProjSet.func);
[Subj, nSubj]=callsubj27;
Respath=ProjSet.Respath;
BrRDM={'Mtv_DMa24'};

%% Turn on the parallel computing
pflag = gcp('nocreate');
if isempty(pflag)
    poolsize=0;
else
    poolsize=pflag.NumWorkers;
end

%% Main
sesIdx=[1:NSes];

for ses=1:NSes
    Cntpart(ses,:)=sesIdx(~ismember(sesIdx,ses)); %other sessions to compare the activity patterns to thoes of the target session
end

for Bri=1:numel(BrRDM)
    disp(BrRDM{Bri});
    savepath=[Respath, BrRDM{Bri}, filesep, 'Vols'];
    if ~exist(savepath,'dir')
        mkdir(savepath);
    end
    subjectPar(Subj, nSubj, BrRDM{Bri}, NSes, Cntpart, savepath);
end %Brain model
delete(gcp('nocreate'));

%% option 1 save RDM
function subjectPar(Subj, nSubj, BrRDM,  NSes, Cntpart, savepath)
    for Sbi=1:nSubj
        disp(Subj{Sbi});
        Lv1path='/mnt/datashare/Pr3.SS/Imaging/Analysis/Analysis Lv1/';
        RDMpath=[Lv1path, BrRDM, filesep, Subj{Sbi}];
        spmTlist=dir([RDMpath, filesep, 'spmT*.nii']);
        spmTname={spmTlist(1:numel(spmTlist)).name};
        Tmat = extractT(RDMpath, spmTname); % num spmT x num voxels 
        RDMs = makeRDM(Tmat, NSes, Cntpart);               
        save([savepath, filesep, Subj{Sbi}, '_sL_RDM.mat'],'RDMs');
    end %Sbj number
end

function Tmat = extractT(RDMpath, spmTname)
    Respath='/mnt/datashare/Pr3.SS/Imaging/Analysis/Analysis RSA/';
	Maskpath=[Respath, 'inclMask.nii']; % Mask: 1 inside the brain, 0 otherwise
    Vmask=spm_read_vols(spm_vol(Maskpath));
    mask=nan(size(Vmask));
    mask(Vmask==1)=1;
    parfor spmTi=1:numel(spmTname)
        spmTpath=[RDMpath, filesep, spmTname{spmTi}];
        data = spm_read_vols(spm_vol(deblank(spmTpath)));
        maskeddata=data.*mask;
	zdata=(maskeddata-nanmean(maskeddata(:))) / nanstd(maskeddata(:));
	%zdata=nanzscore(data.*mask);
	Tmat(spmTi,:) = zdata(:);
    end
end

%% makeRDM option1, if use sL3 to compute tau with many differnt model (save RDM)
function [RDMs] = makeRDM(Tmat, Nses, Cntpart)%[RDMs, shuffleRDM] = makeRDM(Tmat) %if use sL3 to compute tau with many differnt model

    Respath='/mnt/datashare/Pr3.SS/Imaging/Analysis/Analysis RSA/'; %maskpath in RSA_sL0_GenSearchlight
	ROIpath=[Respath, 'Mask29/Mask29_vol/searchlight_100.mat'];
    load(ROIpath); % LI, id of voxels (num voxel x 100 voxels in each searchlight)
    nROIs=size(LI,1);
    
    parfor LIi=1:nROIs 
        SingleSubVol=Tmat(:,LI{LIi});
        nCondition=size(Tmat,1)/Nses;
        RDMs0=cutTmatperses(Tmat,Nses,nCondition,Cntpart, LI{LIi});
        RDMs(LIi) =RDMs0;
    end
end

function [RDMs]=cutTmatperses(Tmat, Nses, nCondition, Cntpart,  LI) % Cutting Tmap per session and Providing two Tmaps to compare (for CrossEuclidean distances)

    for ses=1:Nses                                                     
        SingleSubVol0=Tmat(nCondition*(ses-1)+1:nCondition*ses, LI); % current session {RDMsize x 100} x Nses
        [SingleSubVol{ses}]=ZinROI(nCondition, SingleSubVol0);
    end
    for  ses=1:Nses                                                     
        cntpartses{ses}=CntptVols(SingleSubVol, [Cntpart(ses,:)]); % other session to compare
    end
    for ses=1:Nses                                                     
        RDMeuc0(ses,:,:) =pdist2(SingleSubVol{ses}, cntpartses{ses} , 'euclidean'); %CrossEuclidean
    end
        RDMs.RDMeuc=squeeze(mean(RDMeuc0));
end

function [Cnttemp]=CntptVols(SingleSubVol, Cnti)
    for cnt=1:numel(Cnti)                                                     
        Cnttemp0(cnt,:,:)=SingleSubVol{Cnti(cnt)}; 
    end
    if size(Cnttemp0,1)==1
        Cnttemp=squeeze(Cnttemp0);  % When nSes=2
    else
        Cnttemp=squeeze(mean(Cnttemp0)); % When nSes>2 then average the patterns of other sessions
    end
end

function [SingleSubVol]=ZinROI(nCondition, SingleSubVol0) % zscore
    for i=1:nCondition
        SingleSubVol(i,:)=(SingleSubVol0(i,:) - nanmean(SingleSubVol0(i,:)) ) / nanstd(SingleSubVol0(i,:));
    end
end
