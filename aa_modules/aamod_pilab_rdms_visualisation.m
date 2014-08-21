% Visualise a disvol.
% [aap,resp]=aamod_pilab_rdmvisualisation(aap,task,subj)
function [aap,resp]=aamod_pilab_rdmvisualisation(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get mean data RDM
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
        disvol = loadbetter(vpath);
        % get stimuli
        spath = aas_getfiles_bystream(aap,subj,'pilab_stimuli');
        stimuli = loadbetter(spath);
        
        % prepare output dirs
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);
        figdir = fullfile(resdir,'figures');
        mkdirifneeded(figdir);

        ts = aap.tasklist.currenttask.settings;
        % quick check to avoid accidentally computing this on a full
        % searchlight disvol
        assert(disvol.nfeatures <= ts.maxn,...
            'number of ROIs exceed maxn (%d)',ts.maxn);

        plotrdms_batch('data',disvol.data,'labels',stimuli,'figdir',figdir,...
            'cmap',ts.cmap,'nrows',ts.nrows);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
