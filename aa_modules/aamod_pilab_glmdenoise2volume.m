% Fit with GLMdenoise from KK.

function [aap,resp]=aamod_pilab_glmdenoise2volume(aap,task,subj)

resp='';

switch task
    case 'doit'
        ts = aap.tasklist.currenttask.settings;

        % glmdenoise result (assume vector form rather than 3D vols)
        glpath = aas_getfiles_bystream(aap,subj,'glmdenoise_results');
        % one per split now
        nsplit = size(glpath,1);

        % load mask
        mpath = aas_getfiles_bystream(aap,subj,'pilab_mask');
        V = spm_vol(mpath(1,:));
        mask = spm_read_vols(V) > 0;

        badinds = cell(nsplit,1);
        for sp = 1:nsplit
            fprintf('converting glmdenoise result %d of %d...\n',...
                sp,nsplit);
            results = loadbetter(glpath(sp,:));

            if ts.tmap
                if ts.poolerrors
                    % compute combined error estimates
                    se = repmat(mean(results.modelse{2},2),...
                        [1 size(results.modelse{2},2)]);
                else
                    se = results.modelse{2};
                end
                assert(~all(se(:)==0),'bad standard errors - no bootstrapping?');
                estimates = (results.modelmd{2} ./ se)';
            else
                estimates = results.modelmd{2}';
            end

            % restrict to conditions of interest
            if ~isempty(ts.ignorelabels)
                %assert(issorted(results.regnames),...
                %    'names must be sorted to ignore labels');
                %lastwarn('');
                % experimental untested fix
                [validnames,coninds] = setdiff(results.regnames,...
                    ts.ignorelabels,'stable');
                [b,id] = lastwarn;
                assert(~isequal(id,...
                    'MATLAB:CELL:SETDIFF:RowsFlagIgnored'),...
                    'bad setdiff behavior: old Matlab version?');
                estimates = estimates(coninds,:);
            else
                validnames = results.regnames;
            end

            % restrict to voxels with reasonable r2 and EPI intensity (by same
            % brightness threshold as in GLMdenoise). Note that since
            % glmdenoise almost invariably brings most voxels up to r2 > 0,
            % thresholding by >0 has little effect in practice.
            r2ok = results.R2 > ts.minr2;
            if ts.maskbybright
                % exclude low EPI intensity voxels
                r2ok = r2ok & results.bright;
            end
            if ts.maskbymean
                % exclude voxels where the mean response across conditions
                % is below baseline
                r2ok = r2ok & (mean(estimates,1)>0)';
            end
            if ts.nanmask
                % NaN out bad features - allows different bad features on
                % different splits
                fprintf('NaN masking.\n');
                estimates(:,~r2ok) = NaN;
            end

            % store bad features across splits
            badinds{sp} = ~r2ok';
            nbad = sum(~r2ok);
            fprintf('%d bad features in this split (%.2f%% of total)\n',...
                nbad, 100 * nbad/numel(r2ok));

            % make volume of this split
            tempvol = SPMVolume(estimates,mask,'header',V,'metasamples',...
                struct('labels',{validnames'},'chunks',...
                ones(length(validnames),1)*sp));
            if sp == 1
                vol = tempvol;
            else
                vol = [vol; tempvol];
            end
        end

        % remove hopeless cases - features that are bad across all splits
        badmat = vertcat(badinds{:});
        goodfeatures = ~all(badmat,1);
        vol = vol(:,goodfeatures);

        % update mask
        % indices for mapping glmdenoise results to vol
        maskinds = find(mask);
        masknbefore = length(maskinds);
        mask(maskinds(~goodfeatures)) = false;
        maskinds = find(mask);
        masknafter = length(maskinds);
        nremoved = masknbefore - masknafter;
        if nremoved>0
            fprintf('removed %d bad features (%.2f%% of total)\n',...
              nremoved,100 * nremoved / masknbefore);
        end
        % save updated mask
        V.fname = mpath(1,:);
        spm_write_vol(V,mask);
        aap=aas_desc_outputs(aap,subj,'pilab_mask',mpath);

        outdir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(outdir);
        outpath = fullfile(outdir,'glmdenoisevol.mat');
        save(outpath,'vol','-v7');
        aap=aas_desc_outputs(aap,subj,'pilab_volume',outpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
