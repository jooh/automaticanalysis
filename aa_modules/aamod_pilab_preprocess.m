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

        if ~isempty(ts.setclass)
            fprintf('setting data to %s\n',ts.setclass);
            epivol.data = feval(ts.setclass,epivol.data);
            designvol.data = feval(ts.setclass,designvol.data);
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
            targetlabels,ignorelabels);
        if ~isequal(coninds,1:nlabels)
            covariates = designvol.data(:,ignoreinds);
            designvol = designvol(:,coninds);
            % project out bad cons
            if ~isempty(ignoreinds)
                for c = 1:epivol.desc.samples.nunique.chunks
                    chunkind = epivol.meta.samples.chunks == ...
                        epivol.desc.samples.unique.chunks(c);
                    assert(isequal(chunkind,...
                        designvol.meta.samples.chunks == ...
                        designvol.desc.samples.unique.chunks(c)),...
                        'mismatched chunks in epivol and designvol');
                    epivol.data(chunkind,:) = epivol.data(chunkind,:) - ...
                        covariates(chunkind,:) * ...
                        (covariates(chunkind,:) \ epivol.data(chunkind,:));
                    designvol.data(chunkind,:) = ...
                        designvol.data(chunkind,:) - ...
                        covariates(chunkind,:) * ...
                        (covariates(chunkind,:) \ designvol.data(chunkind,:));
                end
            end
        end

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

        if ~isempty(ts.medianfilter)
            fprintf('median filter (n=%.0f)\n',ts.medianfilter);
            filterbychunk(epivol,'medianfilter',ts.medianfilter);
            filterbychunk(designvol,'medianfilter',ts.medianfilter);
        end

        if ts.sgdetrend
            fprintf('Savitzky-Golay detrend (k=%.0f,f=%.0f)\n',...
                ts.sgolayK,ts.sgolayF);
            sgdetrend(epivol,ts.sgolayK,ts.sgolayF);
            sgdetrend(designvol,ts.sgolayK,ts.sgolayF);
        end

        if ts.zscore
            fprintf('Z-scoring samples\n')
            filterbychunk(epivol,'zscore',[],1);
            filterbychunk(designvol,'zscore',[],1);
        end

        % output
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        out = fullfile(pidir,'epivol.mat');
        save(out,'epivol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_epi',out);
        outdesign = fullfile(pidir,'designvol.mat');
        save(outdesign,'designvol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_design',outdesign);

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sess',...
            outpaths_sessrdms);
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',...
            outpath_mean);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
