function [fmri_spec, condcounter] = Model_Multivariate(subj, nblocks, session, fmri_spec, modelname)
    
    [ProjSet, fs]=CallProjSet;
    data_path = ProjSet.DATApath;
 
    for sess = 1:nblocks

        conditions=modelname(5:end);
        condcounter = 0;

%% Specify regressors
        names=fieldnames( session(sess).(conditions) );  
        for ev=1:numel(fieldnames(session(sess).(conditions) ))
            condcounter = condcounter+1; % condition bases
            FXX=names{ev};
            fmri_spec.sess(sess).cond(1,condcounter).onset             = session(sess).(conditions).(FXX).ons;
            fmri_spec.sess(sess).cond(1,condcounter).tmod              = 0;
            fmri_spec.sess(sess).cond(1,condcounter).pmod              = struct('name', {}, 'param', {}, 'poly', {});
            fmri_spec.sess(sess).cond(1,condcounter).duration          = 2;
            fmri_spec.sess(sess).cond(1,condcounter).orth              = 0; %default in SPM12 is yes which is 1 (orthoginalized)
            %end
            
            for c=1:condcounter
                fmri_spec.sess(sess).cond(c).name = names{c};
            end
        end

%% Add motion regressors
    for sess = 1:nblocks
        fmri_spec.sess(sess).multi           = {''};
        fmri_spec.sess(sess).regress         = struct('name', {}, 'val', {});
        curfd=pwd;
        cd([data_path, subj, '/R', num2str(sess)]);
        rpfile=dir('rp*.txt');
        fmri_spec.sess(sess).multi_reg      = {[data_path, subj, '/R', num2str(sess), fs rpfile.name]}; 
        fmri_spec.sess(sess).hpf            = 128;
        fmri_spec.fact                      = struct('name', {}, 'levels', {});
        fmri_spec.bases.hrf.derivs          = [0 0];
        fmri_spec.global                    = 'None';
        fmri_spec.cvi                       = 'AR(1)';
        cd(curfd);
    end
end
