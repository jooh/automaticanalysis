% export cleaned EPIs and design matrix from glmdenoise to pilab instances.
% [aap,resp]=aamod_pilab_glmdenoisedepi2volume(aap,task,subj)
function [aap,resp]=aamod_pilab_glmdenoisedepi2volume(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='construct glmdenoised EPI and design pilab Volumes';
        
    case 'summary'
        
    case 'report'
        
    case 'doit'
        % glmdenoised results (for design matrix, also get TR from here for
        % frameperoid field in epivol)
        respath = aas_getfiles_bystream(aap,subj,'glmdenoise_results');
        % if separate GLM denoising of different runs
        nsplit = size(respath,1);
        % glmdenoised epis
        epipath = aas_getfiles_bystream(aap,subj,'glmdenoise_epi');
        
        % mask (mainly used so we preserve header in EPI)
        mpath = aas_getfiles_bystream(aap,subj,'pilab_mask');
        V = spm_vol(mpath(1,:));
        mask = spm_read_vols(V) > 0;

        ts = aap.tasklist.currenttask.settings;

        offset = 0;
        badinds = cell(nsplit,1);
        for sp = 1:nsplit
            fprintf('converting glmdenoise EPIs for split %d of %d...\n',...
                sp,nsplit);
            results = loadbetter(respath(sp,:));

            % feature selection
            r2ok = results.R2 > ts.minr2;
            if ts.maskbybright
                % exclude low EPI intensity voxels
                r2ok = r2ok & results.bright;
            end
            if ts.maskbymean
                % exclude voxels where the mean response across conditions
                % is below baseline - note that this is slightly different
                % from the maskbymean option in glmdenoise2volume since we
                % do not first convert to T units or remove ignorelabels.
                % In practice this seems to make little difference.
                r2ok = r2ok & mean(results.modelmd{2},2)>0;
            end

            epi = loadbetter(epipath(sp,:));
            % glmdenoise stores epis with samples in columns and features in
            % rows
            epi = cellfun(@(x)x',epi,'uniformoutput',false);

            if ts.nanmask
                % NaN out bad features - allows different bad features on
                % different splits
                fprintf('NaN masking.\n');
                for ep = 1:length(epi)
                    epi{ep}(:,~r2ok) = NaN;
                end
            end
            % store bad features across splits
            badinds{sp} = ~r2ok';
            nbad = sum(~r2ok);
            fprintf('%d bad features in this split (%.2f%% of total)\n',...
                nbad, 100 * nbad/numel(r2ok));

            nchunks = length(epi);
            % use length to figure out a chunks vector
            chunks = cell2mat(arrayfun(@(x)ones(size(epi{x},1),1)*x,...
                (1:nchunks)','uniformoutput',false)) + offset;
            % keep track of the truly independent splits 
            superchunks = ones(length(chunks),1)*sp;
            % epivol for this split
            tempvol = SPMVolume(epi,mpath(1,:),'frameperiod',...
                results.inputs.tr,'metasamples',struct('chunks',chunks,...
                'superchunks',superchunks));
            % designvol for this split
            design = results.inputs.design;
            if ~strcmp(results.inputs.hrfmodel,'assumeconvolved')
                % need to convolve design matrix - TODO
                error('not yet implemented!')
            end

            extraregs = cell(1,length(design));
            denspec = results.inputs.opt.denoisespec{1};
            assert(length(results.inputs.opt.denoisespec)==1,...
                'only one denoisespec supported');
            if isequal(denspec(1),'0')
                % EPIs have had signal removed
                error('not supported!');
            end
            if isequal(denspec(2),'0')
                % EPIs have had polynomial detrend
                fprintf('projecting out trend polynomials from design\n');
                for r = 1:length(design)
                    nvol = size(design{r},1);
                    extraregs{r} = cat(2,extraregs{r},...
                        constructpolynomialmatrix(nvol,...
                        0:results.inputs.opt.maxpolydeg(r)));
                end
            end
            if isequal(denspec(3),'0')
                % EPIs have had extraregressors removed
                fprintf('projecting out extra regressors from design\n');
                for r = 1:length(design)
                    extraregs{r} = cat(2,extraregs{r},...
                        results.inputs.opt.extraregressors{r});
                end
            else
                assert(isequal(results.inputs.opt.extraregressors{:})&&...
                    isempty(results.inputs.opt.extraregressors{1}),...
                    'EPIs must have extraregressors regressed out');
            end
            if isequal(denspec(4),'0')
                % EPIs have been denoised - design matrix needs this too
                fprintf('projecting out noise PCs from design\n');
                for r = 1:length(design)
                    extraregs{r} = cat(2,extraregs{r},...
                        results.pcregressors{r}(:,1:results.pcnum));
                end
            end
            if isequal(denspec(5),'0')
                % EPIs have had residuals removed
                error('not supported')
            end
            % project if necessary
            if ~isempty(extraregs{1})
                for r = 1:length(design)
                    design{r} = projectionmatrix(extraregs{r}) * design{r};
                end
            end
            tempdes = Volume(vertcat(design{:}),...
                'metafeatures',struct('labels',{results.regnames}),...
                'metasamples',struct('chunks',chunks,'superchunks',...
                superchunks));
            if sp == 1
                epivol = tempvol;
                designvol = tempdes;
            else
                epivol = [epivol; tempvol];
                designvol = [designvol; tempdes];
            end
            % update offset (to make unique chunks across splits)
            offset = max(chunks);
        end
        badmat = vertcat(badinds{:});
        if ts.nanmask
            % save memory by removing hopeless cases - features that are
            % bad across all splits.
            goodfeatures = ~all(badmat,1);
        else
            % remove any feature that fails any test in any split
            % (important to get rid of features that may be unreliable
            % across splits)
            goodfeatures = ~any(badmat,1);
        end
        epivol = epivol(:,goodfeatures);

        % update mask
        % indices for mapping glmdenoise results to vol
        maskinds = find(mask);
        masknbefore = length(maskinds);
        mask(maskinds(~goodfeatures)) = false;
        maskinds = find(mask);
        masknafter = length(maskinds);
        nremoved = masknbefore - masknafter;
        if nremoved>0
            fprintf('removed %d bad features from full volume (%.2f%% of total)\n',...
              nremoved,100 * nremoved / masknbefore);
        end
        % save updated mask
        V.fname = mpath(1,:);
        spm_write_vol(V,mask);
        aap=aas_desc_outputs(aap,subj,'pilab_mask',mpath);
        
        % save epi
        outdir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(outdir);
        outpath_epi = fullfile(outdir,'epivol.mat');
        % very likely too big for older Matlab formats
        save(outpath_epi,'epivol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_epi',outpath_epi);
        outpath_design = fullfile(outdir,'designvol.mat');
        save(outpath_design,'designvol','-v7');
        aap=aas_desc_outputs(aap,subj,'pilab_design',outpath_design);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
