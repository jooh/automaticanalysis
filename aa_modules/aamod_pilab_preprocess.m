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

        % de-trend config
        if strcmp(ts.covariatedeg,'adaptive')
            ts.covariatedeg = vol2covdeg(epivol);
        end

        % first high-pass trend removal
        if ~isempty(ts.covariatedeg)
            fprintf('polynomial detrend (degree=%.0f)\n',ts.covariatedeg);
            filterbychunk(epivol,'polydetrend',ts.covariatedeg);
            filterbychunk(designvol,'polydetrend',ts.covariatedeg);
        end

        if ~isempty(ts.medianfilter) && ts.medianfilter
            fprintf('median filter (n=%.0f)\n',ts.medianfilter);
            filterbychunk(epivol,'medianfilter',ts.medianfilter);
            filterbychunk(designvol,'medianfilter',ts.medianfilter);
        end

        if ts.sgdetrend
            fprintf('Savitzky-Golay detrend (k=%.0f,f=%.0f)\n',...
                ts.sgolayK,ts.sgolayF);
            % insure double
            epivol.data = double(epivol.data);
            designvol.data = double(designvol.data);
            sgdetrend(epivol,ts.sgolayK,ts.sgolayF);
            sgdetrend(designvol,ts.sgolayK,ts.sgolayF);
        end

        if ts.zscore
            fprintf('Z-scoring samples\n')
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
                fprintf('projecting out %d covariates\n',...
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
            fprintf('setting data to %s\n',ts.setclass);
            epivol.data = feval(ts.setclass,epivol.data);
            designvol.data = feval(ts.setclass,designvol.data);
        end

        if ~isempty(ts.resortind)
            if ischar(ts.resortind)
                ts.resortind = feval(ts.resortind,designvol.nfeatures);
            end
            fprintf('resorting regressors in designvol.\n');
            designvol = designvol(:,ts.resortind);
        end

        % output
        fprintf('saving...\n');
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
