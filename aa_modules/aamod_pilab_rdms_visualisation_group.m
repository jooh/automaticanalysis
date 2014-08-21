% Visualise a group disvol.
% [aap,resp]=aamod_pilab_rdmvisualisation_group(aap,task)
function [aap,resp]=aamod_pilab_rdmvisualisation_group(aap,task)

resp='';

switch task
    case 'doit'
        % get RDMs
        meanres = loadbetter(aas_getfiles_bystream(aap,'pilab_rdms_rfx'));
        % get stimuli (NB, we use subject 1's stimuli as an example)
        spath = aas_getfiles_bystream(aap,1,'pilab_stimuli');
        stimuli = loadbetter(spath);
        groupres = loadbetter(aas_getfiles_bystream(aap,'pilab_rdms_group'));
        
        % prepare output dirs
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);
        figdir = fullfile(resdir,'figures');
        mkdirifneeded(figdir);

        ts = aap.tasklist.currenttask.settings;
        % quick check to avoid accidentally computing this on a full
        % searchlight disvol
        assert(size(meanres.mean,2) <= ts.maxn,...
            'number of ROIs exceed maxn (%d)',ts.maxn);

        if ~isempty(ts.pluginpath)
            feval(ts.pluginpath,meanres,stimuli,figdir,ts);
        end

        plotrdms_batch('data',meanres.(ts.target),'roinames',...
            meanres.cols_roi,'labels',stimuli,'figdir',figdir,...
            'cmap',ts.cmap,'nrows',ts.nrows,'gridlines',ts.gridlines,...
            'gridcolor',ts.gridcolor,'ranktransform',ts.ranktransform);

        for roi = 1:numel(meanres.cols_roi)
            roidata = squeeze(groupres.t(:,roi,:));
            if ts.ranktransform==1
                roidata = tiedrank(roidata);
            end
            % nb we omit the labels here because this figure has many
            % panels and gets quite messy
            F(roi) = rdmfig(roidata,meanres.z_subject,[],[],'labels',[],...
                'cmap',ts.cmap,'gridlines',ts.gridlines,'gridcolor',...
                ts.gridcolor);
            set(F(roi),'name',['singles ' meanres.cols_roi{roi}]);
        end
        printbyname(F,figdir);
        close(F);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
