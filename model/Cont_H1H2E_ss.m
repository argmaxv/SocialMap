function [spm] = Cont_H1H2E_ss(cont_names, cntmx, spm)
    
% Assigning weights to generate con files and returning 'spm.stats.con.consess'

    motion_reg=6; %6 motion regressors
    
    for n=1:length(cont_names)
        spm.stats.con.consess{n}.tcon.name =cont_names{n};

        elseif strcmp(cont_names{n}, 'F1H2E')==1                                             %  Contrast 1
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'F1H2E')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;

        elseif strcmp(cont_names{n}, 'F2H2E')==1                                             %  Contrast 2
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'F2H2E')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
        
        elseif strcmp(cont_names{n}, 'F1H2Th')==1                                             %  Contrast 3
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'F1H2Th')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
            
        elseif strcmp(cont_names{n}, 'F2H1Th')==1                                             %  Contrast 4
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'F2H1Th')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;

    % Suppression
        elseif strcmp(cont_names{n}, 'FSNotHB')==1                                             %  Contrast 5
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'FSnotHB')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
        
        elseif strcmp(cont_names{n}, 'FSHB1')==1                                             %  Contrast 6
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'FSisH1')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
        
        elseif strcmp(cont_names{n}, 'FSHB2')==1                                             %  Contrast 7
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'FSisH2')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
        
        elseif strcmp(cont_names{n}, 'FSHB2vsNot')==1                                        %  Contrast 8
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, 'FSisH2')))=-1; % put -1 for regressor of interests because it is for repetition suppression 
            cntw(1,find(strcmp(cntmx, 'FSnotHB')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
        else
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            cntw(1,find(strcmp(cntmx, cont_names{n})))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
        end
        spm.stats.con.consess{n}.tcon.sessrep = 'repl'; %replicates over sessions
    end
    spm.stats.con.delete = 1;
end
