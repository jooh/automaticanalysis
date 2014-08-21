% Import predictor RDMs for RSA.
%
% [aap,resp]=aamod_pilab_importrsapredictors(aap,task,subj)
function [aap,resp]=aamod_pilab_importrsapredictors(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;

        % get stimuli
        spath = aas_getfiles_bystream(aap,subj,'pilab_stimuli');
        stimuli = loadbetter(spath);

        ts = aap.tasklist.currenttask.settings;

        rdms = feval(ts.predictorfun,subname,stimuli);

        % save and describe
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);
        outpath = fullfile(pidir,'pilab_rsapredictors.mat');

        if ~isempty(ts.resortind)
            if ischar(ts.resortind)
                ts.resortind = feval(ts.resortind,numel(stimuli));
            end
            rdmmat = asrdmmat(rdms);
            rdmmat = rdmmat(ts.resortind,ts.resortind,:);
            % assume the stimuli are already sorted appropriately
            for c = 1:length(rdms)
                rdms(c).RDM = rdmmat(:,:,c);
            end
        end
        save(outpath,'rdms');
        aap = aas_desc_outputs(aap,subj,'pilab_rsapredictors',outpath);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
