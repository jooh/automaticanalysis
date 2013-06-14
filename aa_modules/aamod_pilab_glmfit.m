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
        nchunks = epivol.desc.samples.nunique.chunks;

        % detrend before splitting
        if ~isempty(ts.sgolayK)
            assert(isempty(ts.covariatedeg),['probably should not ' ...
                'detrend and include trend covariates at the same time']);
            fprintf('detrending (K=%d,F=%d)\n',...
                ts.sgolayK,...
                ts.sgolayF);
            epivol.sgdetrend(ts.sgolayK,...
                ts.sgolayF);
            designvol.sgdetrend(ts.sgolayK,...
                ts.sgolayF);
        end

        % configure split
        if isempty(ts.split)
            % just one global split
            split = ones(1,nchunks);
        elseif ischar(ts.split)
            split = eval(ts.split);
        end
        usplit = unique(split);
        nsplit = length(usplit);

        % find a reasonable polynomial degree
        if strcmp(ts.covariatedeg,'adaptive')
            covdegs = NaN([1 nchunks]);
            for c = 1:nchunks
                nvol = sum(epivol.meta.samples.chunks == ...
                    epivol.desc.samples.unique.chunks(c));
                % Kay's rule for deciding on n covariates based on run
                % duration.
                covdegs(c) = round(nvol * epivol.frameperiod / 60 / 2);
            end
            assert(~any(isnan(ts.covariatedeg)),...
                'failed to identify covariatedeg');
            % if we found a single deg that works for all runs, life is
            % easy
            if all(covdegs(1)==covdegs)
                ts.covariatedeg = covdegs(1);
                fprintf('adaptively selected covariate degree: %d\n',...
                    ts.covariatedeg);
            else
                % otherwise we need to add new functionality to CovGLM to
                % support different covariatedeg for different runs.
                error(['Adaptive covariate deg selection failed. ' ...
                    'Different covariates selected per run: ' ...
                    mat2str(covdegs)]);
            end
        end

        % find correct labels
        nlabels = length(designvol.desc.features.unique.labels);
        if isempty(ts.targetname)
            coninds = 1:nlabels;
        else
            coninds = findStrInArray(designvol.desc.features.unique.labels,...
                ts.targetname);
        end
        
        if ~isempty(ts.ignorelabels)
            ignoreinds = findStrInArray(designvol.desc.features.unique.labels,...
                ts.ignorelabels);
            coninds = setdiff(coninds,ignoreinds);
        end
       
        assert(~isempty(coninds),...
            'found no labels matching targetname/ignorelabels %s/%s',...
            ts.targetname,ts.ignorelabels);
        ncon = length(coninds);
        fprintf('storing estimates for %d conditions\n',ncon);

        % empty cells don't get read properly
        if isempty(ts.glmvarargs)
            ts.glmvarargs = {};
        end

        % generate separate estimates for each split
        datcell = cell(nsplit,1);
        for s = 1:nsplit
            fprintf('fitting split %d of %d...\n',s,nsplit);
            % split-specific data
            splitind = split==usplit(s);
            splitepi = epivol.selectbymeta('chunks',...
                epivol.desc.samples.unique.chunks(splitind));
            splitdesign = designvol.selectbymeta('chunks',...
                designvol.desc.samples.unique.chunks(splitind));

            % now mask out any NaN features
            nanmask = ~any(isnan(splitepi.data),1);
            % if you haven't NaNed out all volumes for a given feature
            % something is likely broken
            assert(all(all(isnan(splitepi.data(:,~nanmask)))),'inconsistent NaN mask detected');
            nnans = sum(~nanmask);
            if ~all(nanmask)
                fprintf(['removed %d NaN features from analysis ' ...
                    '(%.2f%% of total).\n'],nnans,...
                    100*(nnans/length(nanmask)));
                splitepi = splitepi(:,nanmask);
            end
            % implement some kind of glm
            glm = vol2glm(splitdesign,splitepi,ts.glmclass,...
                ts.covariatedeg,ts.glmvarargs{:});

            % use bootstrap or plain tmapped estimates?
            if ~isempty(ts.nboot) && ts.nboot > 0
                % boot
                [estimates,sterrs] = glm.bootstrapestimate(ts.nboot);
                % optionally convert to variance units
                if ts.tmap
                    estimates = estimates ./ sterrs;
                end
                % restrict to conditions of interest
                estimates = estimates(coninds,:);
            else
                % parametric estimates
                if ts.tmap
                    % compute T contrast for each condition of interest v
                    % baseline
                    labelinds = splitdesign.desc.features.inds.labels;
                    estimates = NaN([ncon glm(1).nfeatures]);
                    regs = zeros(1,glm(1).npredictors);
                    for t = 1:ncon
                        cv = regs;
                        cv(labelinds==coninds(t)) = 1;
                        estimates(t,:) = glm.tmap(cv);
                    end
                else
                    % just get parameter estimate (ie, mean)
                    estimates = glm.fit;
                    estimates = estimates(coninds,:);
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
                struct('chunks',ones(ncon,1)*usplit(s)','labels',...
                {splitdesign.desc.features.unique.labels(coninds)'}));
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
