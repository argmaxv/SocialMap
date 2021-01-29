function [fmri_spec, condcounter] = Model_H1H2E_ss(subj, nblocks, session, fmri_spec)
    
    % The model specification in which the brain activity is modled by the Euclidean distances between F2-Hub1 and F1-Hub2.
    % For testing other models (univariate)
    % 1. Change Level1.m to call differnt Model.
    % 2. Change the name of this function (e.g. Model_H1H2DI_ss)
    % 3. Comment the lines 45-48 and Uncomment the lines associated with the model that you want to test (e.g. lines 39-42 for Model_H1H2DI_ss)
    
    [ProjSet, fs]=CallProjSet;
    data_path = ProjSet.DATApath;

    for sess = 1:nblocks
        
        conditions = {'F0', 'F1', 'F2', 'FSnotHB','FSisH1','FSisH2', 'Btn'};
        cond.conditions = conditions;
        onsetlist=fieldnames(session(sess));
        condcounter = 0;   
        correctTrial=ismember(session(sess).F2.ons, session(sess).F2c.ons);
%% Specify regressors
        for ev=1:length(conditions)
            condcounter = condcounter+1;    

% Parametric regressors
            if strcmp(conditions{ev}, 'F1')==1
                fmri_spec.sess(sess).cond(1,condcounter).onset             = session(sess).(conditions{ev}).ons; 
                fmri_spec.sess(sess).cond(1,condcounter).tmod              = 0;
                %fmri_spec.sess(sess).cond(1,condcounter).pmod                                         = struct('name', {}, 'param', {}, 'poly', {});
                idpm=1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)                 = struct('name', 'F1RC', 'param', session(sess).(conditions{ev}).pm(:,1), 'poly', 1);
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1IC', 'param', session(sess).(conditions{ev}).pm(:,2) , 'poly', 1); 
%                 
            elseif strcmp(conditions{ev}, 'F2')==1
                fmri_spec.sess(sess).cond(1,condcounter).onset             = session(sess).(conditions{ev}).ons(correctTrial); 
                fmri_spec.sess(sess).cond(1,condcounter).tmod              = 0;
                %fmri_spec.sess(sess).cond(1,condcounter).pmod                                         = struct('name', {}, 'param', {}, 'poly', {});
                idpm=1;            fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F2RC', 'param', session(sess).(conditions{ev}).pm(correctTrial,1), 'poly', 1);
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F2IC', 'param', session(sess).(conditions{ev}).pm(correctTrial,2) , 'poly', 1);

                % H1H2DI
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1H2D', 'param', session(sess).(conditions{ev}).pm(correctTrial,7) , 'poly', 1);
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F2H1D', 'param', session(sess).(conditions{ev}).pm(correctTrial,8) , 'poly', 1);
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1H2I', 'param', session(sess).(conditions{ev}).pm(correctTrial,21) , 'poly', 1);
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F2H1I', 'param', session(sess).(conditions{ev}).pm(correctTrial,22) , 'poly', 1);
            
                %H1H2ETh
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1H2E', 'param', session(sess).(conditions{ev}).pm(correctTrial,9) , 'poly', 1);
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F2H1E', 'param', session(sess).(conditions{ev}).pm(correctTrial,10) , 'poly', 1);
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1H2Th', 'param', session(sess).(conditions{ev}).pm(correctTrial,23) , 'poly', 1);
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F2H1Th', 'param', session(sess).(conditions{ev}).pm(correctTrial,24) , 'poly', 1);
                
                % F1F2DI
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1F2D', 'param', session(sess).(conditions{ev}).pm(correctTrial,3) , 'poly', 1);
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F2F1I', 'param', session(sess).(conditions{ev}).pm(correctTrial,13) , 'poly', 1);

                % F1F2ETh
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1F2E', 'param', session(sess).(conditions{ev}).pm(correctTrial,6) , 'poly', 1);
                %idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'F1F2Th', 'param', session(sess).(conditions{ev}).pm(correctTrial25) , 'poly', 1);
                                       
            else 
                fmri_spec.sess(sess).cond(1,condcounter).onset= session(sess).(conditions{ev}).ons; 
                fmri_spec.sess(sess).cond(1,condcounter).tmod = 0;
                fmri_spec.sess(sess).cond(1,condcounter).pmod= struct('name', {}, 'param', {}, 'poly', {});
            end

%% Duration            
            %if isfield(session(sess).(conditions{ev}),'dur')==1
            if strcmp(conditions{ev}, 'Btn')==1 
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 0; %session(sess).(conditions{ev}).dur;
            else
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 2;
            end
            
            fmri_spec.sess(sess).cond(1,condcounter).orth                = 0; %default in SPM12 is yes which is 1 (orthoginalized)
        end
        
        for c=1:condcounter
            fmri_spec.sess(sess).cond(c).name = conditions{c};
        end
        
    end
%clear session;
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
