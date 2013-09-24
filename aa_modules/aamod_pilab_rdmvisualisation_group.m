% Visualise a group disvol.
% [aap,resp]=aamod_pilab_rdmvisualisation_group(aap,task)
function [aap,resp]=aamod_pilab_rdmvisualisation_group(aap,task)

resp='';

switch task
    case 'doit'
        % get mean data RDM
        vpath = aas_getfiles_bystream(aap,'pilab_data_rdms_group_mean');
        disvol = loadbetter(vpath);
        % get stimuli (NB, we use subject 1's stimuli as an example)
        spath = aas_getfiles_bystream(aap,1,'pilab_stimuli');
        stimuli = loadbetter(spath);
        
        % prepare output dirs
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);
        figdir = fullfile(resdir,'figures');
        mkdirifneeded(figdir);

        ts = aap.tasklist.currenttask.settings;
        % quick check to avoid accidentally computing this on a full
        % searchlight disvol
        assert(disvol.nfeatures <= ts.maxn,...
            'number of ROIs exceed maxn (%d)',ts.maxn);

        plotrdms_batch(disvol,stimuli,'figdir',figdir,...
            'cmap',ts.cmap,'nrows',ts.nrows,'alphacolor',ts.alphacolor,...
            'mdstimsize',ts.mdstimsize);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
