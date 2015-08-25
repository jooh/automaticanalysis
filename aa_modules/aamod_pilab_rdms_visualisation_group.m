% Visualise a group disvol.
% [aap,resp]=aamod_pilab_rdmvisualisation_group(aap,task)
function [aap,resp]=aamod_pilab_rdmvisualisation_group(aap,task)

resp='';

switch task
    case 'doit'
        % get RDMs
        meanres = loadbetter(aas_getfiles_bystream(aap,'pilab_result_rfx'));
        % get stimuli (NB, we use subject 1's stimuli as an example)
        spath = aas_getfiles_bystream(aap,1,'pilab_stimuli');
        stimuli = loadbetter(spath);
        ts = aap.tasklist.currenttask.settings;
        rowlabels = stimuli;
        collabels = stimuli;
        if ~isempty(ts.rowind)
            if ischar(ts.rowind)
                ts.rowind = eval(ts.rowind);
            end
            rowlabels = stimuli(ts.rowind);
        end
        if ~isempty(ts.colind)
            if ischar(ts.colind)
                ts.colind = eval(ts.colind);
            end
            collabels = stimuli(ts.colind);
        end

        groupres = loadbetter(aas_getfiles_bystream(aap,'pilab_result_group'));
        
        % prepare output dirs
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        figdir = fullfile(pidir,'figures');
        mkdirifneeded(figdir);

        % quick check to avoid accidentally computing this on a full
        % searchlight disvol
        assert(size(meanres.(ts.mtarget),2) <= ts.maxn,...
            'number of ROIs exceed maxn (%d)',ts.maxn);

        if ~isempty(ts.pluginpath)
            if ~iscell(ts.pluginpath)
                ts.pluginpath = {ts.pluginpath};
            end
            for n = 1:numel(ts.pluginpath)
                feval(ts.pluginpath{n},meanres,groupres,stimuli,figdir,ts);
            end
        end

        if ts.runstandardplots
            plotrdms_batch('data',meanres.(ts.mtarget),'roinames',...
                meanres.cols_roi,'ylabels',rowlabels,'xlabels',collabels,...
                'figdir',figdir,...
                'cmap',ts.cmap,'nrows',ts.nrows,'gridlines',ts.gridlines,...
                'gridcolor',ts.gridcolor,'ranktransform',ts.ranktransform);

            for roi = 1:numel(meanres.cols_roi)
                roidata = squeeze(groupres.d(:,roi,:));
                if ts.ranktransform==1
                    roidata = ranktrans(roidata);
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
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
