% visualise RSA RFX results.
%
% [aap,resp]=aamod_pilab_rsa_visualisation_rfx(aap,task)
function [aap,resp]=aamod_pilab_rsa_visualisation_rfx(aap,task)

resp='';

switch task
    case 'doit'
        % get the results
        meanres = loadbetter(aas_getfiles_bystream(aap,'pilab_result_rfx'));
        groupres = loadbetter(aas_getfiles_bystream(aap,...
            'pilab_result_group'));
        ts = aap.tasklist.currenttask.settings;

        % get RDMs
        nsub = length(aap.acq_details.subjects);
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
            if isstruct(rdms{1})
                meanpredictors(p).RDM = rdms{1};
            else
                meanpredictors(p).RDM = matmean(rdms{:});
            end
        end

        % prepare output
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        figdir = fullfile(pidir,'figures');

        arglist = structfields2varargs(ts.roidataargs);
        if ~isempty(ts.pluginpath)
            feval(ts.pluginpath,figdir,meanres,groupres,meanpredictors,arglist{:});
        end

        % standard plots
        if ts.runstandardplots
            handles = roidata2figure(meanres,groupres,arglist{:});
            printbyname([handles.figure],figdir);
            close([handles.figure]);
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
