function smoothingnii(FWHMmm, prefix, pathfile, mask)
    fprintf('\n Smoothing with a %gmm FWHM Gaussian kernel \n',FWHMmm);
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(pathfile);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHMmm FWHMmm FWHMmm];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = mask;
    matlabbatch{1}.spm.spatial.smooth.prefix = prefix;
    spm_jobman('run',matlabbatch);
end
