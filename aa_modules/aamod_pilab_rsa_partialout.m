% partial out some model RDM from the data + model RDMs (so subsequent fits
% in aamod_pilab_rsa are partial rho/r).
%
% [aap,resp]=aamod_pilab_rsa_partialout(aap,task,subj)
function [aap,resp]=aamod_pilab_rsa(aap,task,subj)

resp='';

switch task
    case 'doit'
        ts = aap.tasklist.currenttask.settings;

        % get data RDMs
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
        disvol = loadbetter(vpath);
        sesspath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_sess');
        splitdisvolcell = arrayfun(@(x)loadbetter(sesspath(x,:)),...
            1:size(sesspath,1),'uniformoutput',0);

        % predictor RDMs
        predictpath = aas_getfiles_bystream(aap,subj,...
            'pilab_rsapredictors');
        if isempty(ts.predictorfun)
            predictors = loadbetter(predictpath);
        else
            logstr('using custom predictor RDMs from %s\n',...
                ts.predictorfun);
            predictors = feval(ts.predictorfun);
        end
        if ~isempty(ts.selectpredictorinds)
            assert(isempty(ts.removepredictorinds),...
                'cannot both select and remove predictor inds');
            if ischar(ts.selectpredictorinds)
                ts.selectpredictorinds = eval(ts.selectpredictorinds);
            end
            predictors = predictors(ts.selectpredictorinds);
        end
        predictors(ts.removepredictorinds) = [];

        if isempty(ts.rsaclassargs) && ~iscell(ts.rsaclassargs)
            ts.rsaclassargs = {};
        end

        assert(~isempty(ts.partialpredictors),...
            'must specify partialpredictors');
        [~,partialind] = intersect({predictors.name},...
            ts.partialpredictors);
        assert(~isempty(partialind),'no match for partialpredictors');
        partialrdms = predictors(partialind);
        % doesn't make sense to partial the partial predictors 
        predictors(partialind) = [];

        % we want to drop any predictors who have NaNs that do
        % not match partialpredictors.
        nanmask = isnan(asrdmmat(partialrdms));
        % check that the nans are consistent among the partialrdms
        ok = bsxfun(@eq,nanmask(:,:,1),nanmask);
        assert(all(ok(:)),'mismatched nans across partialpredictors');
        % then drop any non-partial predictors with the same issue
        oktargets = arrayfun(@(x)isequal(nanmask(:,:,1),...
            isnan(asrdmmat(predictors(:,x)))),1:numel(predictors));
        if ~all(oktargets)
            logstr(['removing %d/%d predictors with mismatched nans ' ...
                'compared to partialpredictors\n'],sum(~oktargets(:)),...
                numel(oktargets));
            predictors(~oktargets) = [];
        end

        % get the residual predictors
        model_predict = feval(ts.rsaclass,partialrdms,predictors,...
            ts.rsaclassargs{:});
        res_predict = asrdmmat(residualsfull(model_predict));
        predictors = rdm2struct(res_predict,{predictors.name});
        save(predictpath,'predictors');

        % residual disvol
        model_disvol = feval(ts.rsaclass,partialrdms,disvol.data);
        disvol.data = residualsfull(model_disvol);
        save(vpath,'disvol');

        % process the splitdisvolcell
        if numel(splitdisvolcell)==1
            sessdisvol = disvol;
            save(sesspath,'sessdisvol');
        else
            for s = 1:numel(splitdisvolcell)
                sessdisvol = splitdisvolcell{s};
                model_sessdisvol = feval(ts.rsaclass,partialrdms,...
                    sessdisvol.data);
                sessdisvol.data = residualsfull(model_sessdisvol);
                save(sesspath(s,:),'sessdisvol');
            end
        end

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sess',sesspath);
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',vpath);
        aap=aas_desc_outputs(aap,subj,'pilab_rsapredictors',predictpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
