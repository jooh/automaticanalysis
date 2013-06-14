% Fit with GLMdenoise from KK.

function [aap,resp]=aamod_firstlevel_glmdenoise(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='run GLMDenoise';
        
    case 'summary'
        
    case 'report'
        
    case 'doit'
        % model (convolved design matrix, so after aamod_firstlevel_design)
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        SPM = loadbetter(spmpath);
        frameperiod = SPM.xY.RT;
        % get mask - memory problems otherwise
        mpath = aas_getfiles_bystream(aap,subj,'epiBETmask');
        % first mask is the brain mask
        V = spm_vol(mpath(1,:));
        mask = spm_read_vols(V) > 0;

        % pass off to new independent function
        ts = aap.tasklist.currenttask.settings;
        % converts and gets rid of NaN / 0 voxels.
        [epi,design,dur,mask,names,extras] = spm2glmdenoise(SPM,mask,...
          ts.ignorelabels,strcmp(ts.hrfmodel,'assumeconvolved'),...
            ts.includeextras);
        % XML parsing probably interprets this as a number
        if ~isempty(ts.opt.denoisespec) && ~ischar(ts.opt.denoisespec)
            ts.opt.denoisespec = sprintf('%05d',ts.opt.denoisespec);
        end

        % optional low pass filter
        if ~isempty(ts.K)
            fprintf('low pass filtering data\n')
            for r = 1:length(epi)
                % NB in GLMdenoise format epi time is in columns rather than in
                % rows
                epi{r} = single(sgolayfilt(double(epi{r}),ts.K,ts.F,[],2));
                % BUT design matrix still has time in rows...
                design{r} = single(sgolayfilt(double(design{r}),ts.K,ts.F,[],1));
            end
        end

        % configure split
        if isempty(ts.split)
            % just one global split
            split = ones(1,length(epi));
        elseif ischar(ts.split)
            split = eval(ts.split);
        end
        usplit = unique(split);
        nsplit = length(usplit);

        % Kendrick's empirical HRFs look very noisy for short stimuli (<3
        % s) and don't work at all for duration 0 (impulse response). In
        % any case they are extremely similar to the spm_hrf (reassuring!).
        % So we are going to use the SPM HRF when durations are 0. The
        % situation is more complicated for longer durations since you'd
        % then need to convolve the spm_hrf to make a predicted HRF for
        % longer durations.
        % (note that the duration input is basically redundant - Kendrick's
        % code only uses the duration for creating the HRF so when this is
        % given (e.g. when 'optimize' and 'hrfknobs' gives an HRF),
        % duration is not used.
        if dur==0 && isempty(ts.hrfknobs) && ~strcmp(ts.hrfmodel,'assumeconvolved')
            fprintf('event duration 0 so using spm HRF\n');
            ts.hrfknobs = normalizemax(spm_hrf(frameperiod));
        end

        if ~isfield(ts.opt,'brainmask') || isempty(ts.opt.brainmask)
            ts.opt.brainmask = mask;
        end

        subdir = aas_getsubjpath(aap,subj);
        outdir = fullfile(subdir,'glmdenoise');
        mkdirifneeded(outdir);

        % fields to write out as diagnostic niftis
        diagnostics = struct('field',{'R2','SNR','noisepool','bright'},...
            'paths',repmat({cell(nsplit,1)},[1 4]));
        ndia = length(diagnostics);

        % run glmdenoise separately on each split
        outpath_models = cell(nsplit,1);
        outpath_results = cell(nsplit,1);
        outpath_epi = cell(nsplit,1);
        for sp = 1:nsplit
            fprintf('glmdenoise split %d of %d...\n',sp,nsplit);
            sessoutdir = fullfile(outdir,sprintf('split%02d',sp));
            mkdirifneeded(sessoutdir);
            figdir = fullfile(sessoutdir,'diagnostic_figures');
            mkdirifneeded(figdir);
            sessinds = find(split==usplit(sp));
            if ts.includeextras
                ts.opt.extraregressors = extras(sessinds);
            end
            [results,denoisedepi] = GLMdenoisedata(design(sessinds),...
                epi(sessinds),dur,frameperiod,ts.hrfmodel,ts.hrfknobs,...
                ts.opt,figdir);
            % add names field
            results.regnames = names;
            % recreate 'bright' logical mask marking voxel exceeding intensity
            % threshold
            results.bright = results.meanvol(:) > (prctile(results.meanvol,...
              results.inputs.opt.brainthresh(1)) * ...
                results.inputs.opt.brainthresh(2));
            % split off very large and fairly irrelevant field in results to
            % prevent ridiculously big MAT files
            models = results.models;
            results.models = [];
            outpath_models{sp} = fullfile(sessoutdir,'results_models.mat');
            save(outpath_models{sp},'models','-v7.3');
            % results struct
            outpath_results{sp} = fullfile(sessoutdir,'results.mat');
            % need 7.3 flag since this may be >2GB
            save(outpath_results{sp},'results','-v7.3');
            % denoised EPIs
            outpath_epi{sp} = fullfile(sessoutdir,'denoisedepi.mat');
            save(outpath_epi{sp},'denoisedepi','-v7.3');
            % write out diagnostic volumes
            for dia = 1:ndia
                diastr = diagnostics(dia).field;
                temp = double(mask);
                temp(mask) = results.(diastr);
                [diagnostics(dia).paths{sp},V.fname] = deal(...
                    fullfile(sessoutdir,[diastr '.nii']));
                V.dt = [spm_type('float32') spm_platform('bigend')];
                spm_write_vol(V,temp);
            end
        end

        % describe outputs - one per split
        aap=aas_desc_outputs(aap,subj,'glmdenoise_results_models',outpath_models);
        aap=aas_desc_outputs(aap,subj,'glmdenoise_results',outpath_results);
        aap=aas_desc_outputs(aap,subj,'glmdenoise_epi',outpath_epi);

        % describe diagnostics
        for dia = 1:ndia
            aap=aas_desc_outputs(aap,subj,['glmdenoise_diagnostic_' ...
                lower(diagnostics(dia).field)],diagnostics(dia).paths);
        end

        % update global mask
        V.fname = mpath(1,:);
        spm_write_vol(V,mask);
        aap=aas_desc_outputs(aap,subj,'epiBETmask',mpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
