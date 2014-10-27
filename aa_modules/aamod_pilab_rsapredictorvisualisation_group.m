% Visualise group predictors.
%
% [aap,resp]=aamod_pilab_rsapredictorvisualisation_group(aap,task)
function [aap,resp]=aamod_pilab_rsapredictorvisualisation_group(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);
        % save results in main module directory
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        % get RDMs
        for s = 1:nsub
            predictors{s} = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsapredictors'));
        end

        % average
        npredict = numel(predictors{1});
        ptest = cellfun(@numel,predictors,'uniformoutput',false);
        assert(isequal(npredict,ptest{:}),...
            'different n predictors across subjects');
        for p = 1:npredict
            meanpredictors(p).name = predictors{1}(p).name;
            ntest = cellfun(@(thisp)thisp(p).name,predictors,...
                'uniformoutput',false);
            assert(isequal(meanpredictors(p).name,ntest{:}),...
                'badly sorted predictors');
            rdms = cellfun(@(thisp)thisp(p).RDM,predictors,...
                'uniformoutput',false);
            meanpredictors(p).RDM = matmean(rdms{:});
        end

        % get stimuli (NB, we use subject 1's stimuli as an example)
        spath = aas_getfiles_bystream(aap,1,'pilab_stimuli');
        stimuli = loadbetter(spath);
        
        % prepare output dirs
        figdir = fullfile(pidir,'figures');
        mkdirifneeded(figdir);

        ts = aap.tasklist.currenttask.settings;

        if ~isempty(ts.pluginpath)
            if ~iscell(ts.pluginpath)
                ts.pluginpath = {ts.pluginpath};
            end
            for n = 1:numel(ts.pluginpath)
                feval(ts.pluginpath{n},meanpredictors,stimuli,figdir,ts);
            end
        end

        plotrdms_batch('data',meanpredictors,'labels',stimuli,...
            'figdir',figdir,...
            'cmap',ts.cmap,'nrows',ts.nrows,'gridlines',ts.gridlines,...
            'gridcolor',ts.gridcolor,'ranktransform',ts.ranktransform);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
