% AA module - add dummy regressors for spike volumes from tsdiffana to
% firstlevel model
% [aap,resp]=aamod_firstlevel_spikes(aap,task,subj)
function [aap,resp]=aamod_firstlevel_spikes(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        %% Add spikes to model
        sessnuminspm=0;
        for sess = aap.acq_details.selected_sessions
            sessnuminspm=sessnuminspm+1;
            % find spikes
            spikepath = aas_getfiles_bystream(aap,subj,sess,'spikes');
            spikes = load(spikepath);
            spikes = spikes.spikes;
            nspikes = length(spikes);
            % make dummy regressors
            zerov = zeros(SPM.nscan(sessnuminspm),1);
            for s = 1:nspikes
                zc = zerov;
                zc(spikes(s)) = 1;
                name = sprintf('spike_%02d',s);
                if ~isfield(SPM,'Sess') || (length(SPM.Sess) < sessnuminspm)
                    SPM.Sess(sessnuminspm) = struct('C',struct('C',[],...
                        'name',[]));
                end
                SPM.Sess(sessnuminspm).C.C = ...
                    [SPM.Sess(sessnuminspm).C.C zc];
                SPM.Sess(sessnuminspm).C.name = ...
                    [SPM.Sess(sessnuminspm).C.name {name}];
            end
            fprintf('added %d spike regressors to session %d\n',nspikes,...
                sess);
        end
        save(spmpath,'SPM');
        % Describe outputs
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',spmpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
