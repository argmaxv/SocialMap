%% 1st level analysis with generating contrasts
% by SPARK 12/28/2017

clear all;
close all;
clc;

%% Setting
[ProjSet, fs, nblocks]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.func);
[subj,nsubj]=callsubj27;

base_path=ProjSet.basepath;
data_path = ProjSet.DATApath;
an_path =ProjSet.ANA1path;
ons_dir=ProjSet.ONSETpath;
addpath ([base_path fs 'Batch' fs]);
addpath ([base_path fs 'Batch' fs 'Model' fs]);

modelspec=1; %0 if you want to generate new contrasts based on what has been modeled before, otherwise 1
contrstgen=1; %0 if you don't need to make con files from the model otherwise 1

design_name = 'H1H2E_ss';
cont_names = {'F1H2E', 'F2H1E', 'F1H2Th', 'F2H1Th',  'FSNotHB', 'FSHB1', 'FSHB2', 'FSHB2vsNot'};

%% Basic model specification
for s = 1:length(subj)
    
    if modelspec==1
                disp(['%%  Starting model ',design_name,' for subject : ', subj{s},'   %%']);
        if ~exist([an_path, design_name, fs, subj{s}, fs],'dir')
            mkdir([an_path, design_name, fs, subj{s}, fs]) 
        end

        cd([an_path,design_name,fs,subj{s},fs]);

        fmri_spec.dir             = {[an_path,design_name,fs,subj{s},fs]}; % directory for SPM.mat file
        fmri_spec.timing.units    = 'secs';
        fmri_spec.timing.RT       = 1.2;        % TR
        fmri_spec.timing.fmri_t   = 38;        % microtime resolution (time-bins per scan) = slices
        fmri_spec.timing.fmri_t0  = 19;       % reference slice from pre-processing
        fmri_spec.mask            = {'/usr/local/MATLAB/spm12/tpm/mask_ICV.nii'}; %selecte mask
        fmri_spec.mthresh         = -Inf; %level of threshold of the mask 
        fmri_spec.volt            = 1;

    % number of scans and data
        for sess = 1:nblocks
            epiimgs =spm_select('List', [data_path,subj{s}, fs, 'R',num2str(sess)], '^swua.*.img');  
            epiimgs=strcat([data_path,subj{s}, fs, 'R',num2str(sess),'/'], epiimgs);
            for epi=1:length(epiimgs)
                fmri_spec.sess(sess).scans{epi, 1}=epiimgs(epi, :);     
            end
        end

    %% Loading mat file extracted from behavior results (it has the matrix called session(sess))
        load ([ons_dir, 'fMRI_Bhv_mat', fs, 'Onset_', subj{s}, '.mat']);
        session=Session.se;

    %% Calling Model specification function
       eval(['[fmri_spec, conditions] = Model_' design_name '(subj{s}, nblocks, session, fmri_spec);']);

    %% Run design specification
        matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
        disp(['%%%%%%%%% Starting design spec for subj: ', subj{s},'%%%%%%%%%']);
        spm_jobman('run',matlabbatch);
        batchname=strcat(fmri_spec.dir{1},design_name,'_model.mat');
        save(batchname,'fmri_spec');
        clear matlabbatch; clear fmri_spec; 

        %% Estimate
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[an_path,design_name,'/',subj{s},'/SPM.mat']};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

        disp(['%%%%%%%%% Starting estimation for subj: ', subj{s},' %%%%%%%%%']);
        spm_jobman('run',matlabbatch);
    end
    
%% Contrast
    if contrstgen==1
    %% Contrast generation

    %Make a folder for the contrast results if not made
        if ~exist([an_path, design_name, fs, subj{s}, fs],'dir')
            mkdir([an_path, design_name, fs, subj{s}, fs])
        end

    %change to new directory for saving
        cd([an_path, design_name, fs, subj{s}, fs])

    %load in lvl1 mat file
        spm.stats.con.spmmat = {[an_path design_name fs subj{s} fs 'SPM.mat']};
        load([design_name,'_model.mat']);
        cntm={}; cntmx={}; im=1; imx=1;
        for cntidx=1:length(fmri_spec.sess(1).cond) %if contrast of session 2 is different from that of session 2 then reprogram the following for loop.

            cntm{im}=fmri_spec.sess(1).cond(cntidx).name;
            im=im+1; %non-parametric regressors

            cntmx{imx}=fmri_spec.sess(1).cond(cntidx).name;
            imx=imx+1; 

            if length(fmri_spec.sess(1).cond(cntidx).pmod)>0
                for cntidy=1:length(fmri_spec.sess(1).cond(cntidx).pmod)
                    cntmx{imx}=fmri_spec.sess(1).cond(cntidx).pmod(cntidy).name; %non + parametric regressors
                    imx=imx+1; 
                end
            else
            end
        end
    %% Contrasts 
    % Calling contrast generate function
        eval(['[spm] = Cont_' design_name '(cont_names, cntmx, spm);']);

    %% Run Contrasts generate Batch
        disp(['%%%%%%%%% Running Contrasts for subj: ', subj{s},'%%%%%%%%%']);
        matlabbatch{1}.spm=spm;
        spm_jobman('run',matlabbatch);
        clear matlabbatch spm Session session epi epiimags fmri_spec
    end % contrastsgen
end%subject loop

% save GLM specification
cd([an_path, design_name, fs]);
glmname=['GLM_' design_name '.mat'];
if exist(glmname) ~= 2
    glmspec.npm_reg=cntm; %non-parametric regressors
    glmspec.pm_reg=cntmx; %Regressors (non+parametric regressors)
    glmspec.cont=cont_names; %Regressor of interests (number is coresponding to the con files)
save(glmname, 'glmspec');
end
