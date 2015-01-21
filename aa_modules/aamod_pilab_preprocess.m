% Preprocess an epivol / designvol pair.
%
% [aap,resp]=aamod_pilab_preprocess(aap,task,subj)
function [aap,resp]=aamod_pilab_preprocess(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        designpath = aas_getfiles_bystream(aap,subj,'pilab_design');
        designvol = loadbetter(designpath);
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);
        ts = aap.tasklist.currenttask.settings;

        if ts.matchn
            nperchunk = arrayfun(@(c)sum(epivol.meta.samples.chunks==c),...
                epivol.desc.samples.unique.chunks);
            targetn = min(nperchunk);
            if all(targetn == nperchunk);
                logstr('all chunks have same number of samples\n');
            else
                logstr('matching nsamples to smallest chunk (%d)\n',...
                    targetn);
                % this might actually be a tad hairy. need to apply the
                % same to design and epi obviously
                goodsamp = false(epivol.nsamples,1);
                for c = 1:epivol.desc.samples.nunique.chunks
                    chunkind = find(epivol.meta.samples.chunks == ...
                        epivol.desc.samples.unique.chunks(c));
                    goodsamp(chunkind(1:targetn)) = true;
                end
                epivol = epivol(goodsamp,:);
                designvol = designvol(goodsamp,:);
                logstr('removed %d samples (%2.0f%% of total)\n',...
                    sum(~goodsamp),100*sum(~goodsamp) / numel(goodsamp));
            end
        end

        % de-trend config
        if strcmp(ts.covariatedeg,'adaptive')
            ts.covariatedeg = vol2covdeg(epivol);
        end

        % first high-pass trend removal
        if ~isempty(ts.covariatedeg)
            logstr('polynomial detrend (degree=%.0f)\n',ts.covariatedeg);
            filterbychunk(epivol,'polydetrend',ts.covariatedeg);
            filterbychunk(designvol,'polydetrend',ts.covariatedeg);
        end

        if ~isempty(ts.medianfilter) && ts.medianfilter
            logstr('median filter (n=%.0f)\n',ts.medianfilter);
            filterbychunk(epivol,'medianfilter',ts.medianfilter);
            filterbychunk(designvol,'medianfilter',ts.medianfilter);
        end

        if ts.sgdetrend
            logstr('Savitzky-Golay detrend (k=%.0f,f=%.0f)\n',...
                ts.sgolayK,ts.sgolayF);
            % insure double
            epivol.data = double(epivol.data);
            designvol.data = double(designvol.data);
            sgdetrend(epivol,ts.sgolayK,ts.sgolayF);
            sgdetrend(designvol,ts.sgolayK,ts.sgolayF);
        end

        if ts.zscore
            logstr('Z-scoring samples\n')
            filterbychunk(epivol,'zscore',[],1);
            filterbychunk(designvol,'zscore',[],1);
        end

        % maybe deal with covariates
        nlabels = length(designvol.desc.features.unique.labels);
        if isempty(ts.targetlabels)
            coninds = 1:nlabels;
        else
            coninds = find(strcmp(designvol.desc.features.unique.labels,...
                ts.targetlabels));
        end
        if ~isempty(ts.ignorelabels)
            ignoreinds = find(strcmp(designvol.desc.features.unique.labels,...
                ts.ignorelabels));
            coninds = setdiff(coninds,ignoreinds,'stable');
        end
        assert(~isempty(coninds),...
            'found no labels matching targetlabels/ignorelabels %s/%s',...
            ts.targetlabels,ts.ignorelabels);
        if ~isequal(coninds,1:nlabels)
            covariates = designvol.data(:,ignoreinds);
            designvol = designvol(:,coninds);
            % project out bad cons
            if ~isempty(ignoreinds)
                logstr('projecting out %d covariates\n',...
                    size(covariates,2));
                for c = 1:epivol.desc.samples.nunique.chunks
                    chunkind = epivol.meta.samples.chunks == ...
                        epivol.desc.samples.unique.chunks(c);
                    assert(isequal(chunkind,...
                        designvol.meta.samples.chunks == ...
                        designvol.desc.samples.unique.chunks(c)),...
                        'mismatched chunks in epivol and designvol');
                    epivol.data(chunkind,:) = projectout(...
                        epivol.data(chunkind,:),covariates(chunkind,:));
                    designvol.data(chunkind,:) = projectout(...
                        designvol.data(chunkind,:),covariates(chunkind,:));
                end
            end
        end

        % now set class last - so maximal precision for pre-processing
        if ~isempty(ts.setclass)
            logstr('setting data to %s\n',ts.setclass);
            epivol.data = feval(ts.setclass,epivol.data);
            designvol.data = feval(ts.setclass,designvol.data);
        end

        if ~isempty(ts.resortind)
            if ischar(ts.resortind)
                ts.resortind = feval(ts.resortind,designvol.nfeatures);
            end
            logstr('resorting regressors in designvol.\n');
            designvol = designvol(:,ts.resortind);
        end

        % output
        logstr('saving...\n');
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        out = fullfile(pidir,'epivol.mat');
        save(out,'epivol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_epi',out);
        outdesign = fullfile(pidir,'designvol.mat');
        save(outdesign,'designvol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_design',outdesign);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
