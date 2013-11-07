% generate a pilab_volume of response estimates from an epivol and
% designvol.
% [aap,resp]=aamod_pilab_glmfit(aap,task,subj)
function [aap,resp]=aamod_pilab_glmfit(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        designpath = aas_getfiles_bystream(aap,subj,'pilab_design');
        designvol = loadbetter(designpath);
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);
        ts = aap.tasklist.currenttask.settings;

        % NB the semantics of this new organisation is quite different. I
        % think we used to allow regressors of no interest in the GLM. Now
        % these are removed before the GLM instance is created.
        predictornames = designvol.meta.features.labels;
        [epicell,designcell] = splitvol(ts.split,epivol,designvol);
        nsplit = length(epicell);
        ncon = designvol.nfeatures;

        if isempty(ts.glmvarargs)
            ts.glmvarargs = {};
        end
        if ~iscell(ts.glmvarargs)
            ts.glmvarargs = {ts.glmvarargs};
        end

        % generate separate estimates for each split
        datcell = cell(nsplit,1);
        for s = 1:nsplit
            nanmask = ~any(isnan(epicell{s}.data),1);
            if ~all(nanmask)
                epicell{s} = epicell{s}(:,nanmask);
            end
            glm = vol2glm(designcell{s},epicell{s},ts.glmclass,...
                ts.glmvarargs{:});

            % use bootstrap or plain tmapped estimates?
            if ~isempty(ts.nboot) && ts.nboot > 0
                % boot
                [estimates,sterrs] = glm.bootstrapestimate(ts.nboot);
                % optionally convert to variance units
                if ts.tmap
                    estimates = estimates ./ sterrs;
                end
                % restrict to conditions of interest
                %estimates = estimates(coninds,:);
            else
                % parametric estimates
                if ts.tmap
                    % compute T contrast for each condition of interest v
                    % baseline
                    estimates = NaN([ncon glm(1).nfeatures]);
                    cmat = eye(ncon);
                    estimates = glm.tmap(cmat);
                else
                    % just get parameter estimate (ie, mean)
                    estimates = glm.fit;
                    %estimates = estimates(coninds,:);
                end
            end

            % re-introduce the NaN'ed out feature. This is necessary to
            % avoid having different number of features and masks for
            % different runs. NaN features then get removed again in e.g.
            % roidata2rdmvol.
            if ~all(nanmask)
                nanmat = repmat(nanmask,[ncon 1]);
                estemp = NaN(size(nanmat));
                estemp(nanmat) = estimates;
                estemp(~nanmat) = NaN;
                estimates = estemp;
            end

            % construct a volume with estimates for this split. Now, the
            % chunk becomes the split (whereas before it was probably the
            % run or sub-run)
            datcell{s} = MriVolume(estimates,epivol,'metasamples',...
                struct('chunks',ones(ncon,1)*s,'labels',...
                {predictornames'}));
        end

        % finally, make a big vol
        vol = vertcat(datcell{:});
        % save and describe
        outdir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpath = fullfile(outdir,'tvol.mat');
        % very likely too big for older Matlab formats
        save(outpath,'vol','-v7');
        aap=aas_desc_outputs(aap,subj,'pilab_volume',outpath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
