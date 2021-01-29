%% 2st level analysis - one sample t test
% by SPARK 12/28/2017

clear all;
close all;
clc

%% Setting
[ProjSet, fs, nblocks]=CallProjSet;
addpath(ProjSet.spmdir);
addpath(ProjSet.func);
[subj,nsubj]=callsubj27;
lvl1folder=ProjSet.ANA1path; 
lvl2folder=ProjSet.ANA2path; 

copyConds = 1;
design_name = 'H1H2E_ss';
cd([lvl1folder, design_name, fs]);
glmlv1=dir('GLM*.mat');
load (glmlv1.name);

for icnt=1:length(glmspec.cont)
    if ~exist([lvl2folder, design_name, fs, glmspec.cont{icnt}, fs, 'con', fs],'dir')
            mkdir([lvl2folder, design_name, fs, glmspec.cont{icnt}, fs, 'con', fs])
    end
end
 
 %copy files over to the level 2 folder
 if copyConds ==1
     %make directories for all subjects
     for s=1:length(subj)
         for icnt=1:length(glmspec.cont)
            confile=spm_select('List', [lvl1folder,design_name, '/', subj{s} ], ['con.*' sprintf('%04d',icnt) '.nii']); % get the name of the original con file 
            confile=cellstr([repmat([lvl1folder,design_name, '/', subj{s},'/'],size(confile,1),1) confile]); % add path to the confile
            copyfile(confile{:}, [lvl2folder, design_name, fs, glmspec.cont{icnt}, fs 'con', fs, subj{s}, '.nii'], 'f'); % copy the confile to the Lv2 folder renaming with subjectID
         end
         
    end
 end   
 
 %pull contrast files
 for icnt=1:length(glmspec.cont)
    lv2dir=[lvl2folder, design_name, fs, glmspec.cont{icnt}, fs];
    factorial_design.dir = {lv2dir};

    for ss=1:length(subj)
      scans{ss}=[lv2dir, fs, 'con' fs subj{ss} '.nii'];
    end
    factorial_design.des.t1.scans=cellstr(scans)';
    
    factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    factorial_design.masking.tm.tm_none = 1;
    factorial_design.masking.im = 1;
    factorial_design.masking.em = {''};
    factorial_design.globalc.g_omit = 1;
    factorial_design.globalm.gmsca.gmsca_no = 1;
    factorial_design.globalm.glonorm = 1;
    disp([num2str(icnt), '.  ', glmspec.cont{icnt}]);
    matlabbatch{1}.spm.stats.factorial_design=factorial_design;    
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {[lv2dir, 'SPM.mat']};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    con.spmmat = {[lv2dir, 'SPM.mat']};
    tconname={[glmspec.cont{icnt} '_pos'], [glmspec.cont{icnt} '_neg']};
    con.consess{1}.tcon.name = tconname{1}; %+'
    con.consess{1}.tcon.weights = 1;
    con.consess{1}.tcon.sessrep = 'none';
    con.consess{2}.tcon.name = tconname{2}; %'-'
    con.consess{2}.tcon.weights = -1;
    con.consess{2}.tcon.sessrep = 'none';
    con.delete = 0;
    matlabbatch{3}.spm.stats.con=con;
    
    results.spmmat ={[lv2dir, 'SPM.mat']};
    for tconN=1:numel(tconname)      
        results.conspec(tconN).titlestr = '';
        results.conspec(tconN).contrasts = tconN;
        results.conspec(tconN).threshdesc = 'none';
        results.conspec(tconN).thresh = 0.001;
        results.conspec(tconN).extent = 0;
        results.conspec(tconN).conjunction = 1;
        results.conspec(tconN).mask.none = 1;
    end
    results.units = 1;
    results.export{1}.ps = true;
    matlabbatch{4}.spm.stats.results=results;
    
    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
 end
 
 
