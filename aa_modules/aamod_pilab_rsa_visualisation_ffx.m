% visualise RSA FFX results.
%
% [aap,resp]=aamod_pilab_rsa_visualisation_ffx(aap,task)
function [aap,resp]=aamod_pilab_rsa_visualisation_ffx(aap,task)

resp='';

switch task
    case 'doit'
        % get the results
        meanres = loadbetter(aas_getfiles_bystream(aap,'pilab_res_ffx'));
        ts = aap.tasklist.currenttask.settings;
        nsub = length(aap.acq_details.subjects);

        % get model RDMs
        for s = 1:nsub
            predictors{s} = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsapredictors'));
        end
        % average
        npredict = numel(predictors{1});
        ptest = cellfun(@numel,predictors,'uniformoutput',false);
        assert(isequal(npredict,ptest{:}),...
            'different n predictors across subjects');
        meanpredictors = predictors{1};
        for p = 1:npredict
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

        arglist = {figdir,meanres,[],meanpredictors,...
            'mtarget',ts.mtarget,...
            'errtarget',ts.errtarget,'ptarget',ts.ptarget,...
            'mlabel',ts.mlabel,'errlabel',ts.errlabel,...
            'pthresh',ts.pthresh};
        if ~isempty(ts.pluginpath)
            % call on mean
            feval(ts.pluginpath,arglist{:});
        end

        if ts.runstandardplots
            % standard plots
            plot_roidata(arglist{:});
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
