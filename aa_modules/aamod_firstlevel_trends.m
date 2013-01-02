% AA module - add trend covariates (Fourier set) to firstlevel model. NB
% this makes high pass filtering redundant.
% [aap,resp]=aamod_firstlevel_trends(aap,task,subj)
function [aap,resp]=aamod_firstlevel_trends(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        sessnuminspm=0;
        for sess = 1:length(SPM.Sess);
            nvol = SPM.nscan(sess);
            ncycles = aap.tasklist.currenttask.settings.ncycles;
            volrange = (1:nvol)';
            doramp = aap.tasklist.currenttask.settings.linearterm;
            % make 2*cycle covariates (sine and cosine) 
            nreg = ncycles*2 + doramp;
            tm = NaN([nvol nreg]);
            names = cell(1,nreg);
            for n = 1:2:(ncycles*2)
                period = nvol / n;
                tm(:,n:n+1) = [sin(volrange/period*2*pi) ...
                    cos(volrange/period*2*pi)];
                names{n} = sprintf('sin%02d',n);
                names{n+1} = sprintf('cos%02d',n);
            end
            % and maybe also a linear ramp
            if aap.tasklist.currenttask.settings.linearterm
                ramp = volrange / nvol;
                ramp = ramp - mean(ramp);
                tm(:,end) = ramp;
                names{end} = 'linear';
            end
            % add to SPM
            SPM.Sess(sess).C.C    = [SPM.Sess(sess).C.C ...
                tm];
            SPM.Sess(sess).C.name = [...
                SPM.Sess(sess).C.name names];
        end
        save(spmpath,'SPM');
        %% Describe outputs
        % Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',spmpath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
