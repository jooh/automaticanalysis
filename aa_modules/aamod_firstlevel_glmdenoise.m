% Fit with GLMdenoise from KK.

function [aap,resp]=aamod_firstlevel_glmdenoise(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='construct EPI and design pilab Volumes';
        
    case 'summary'
        
    case 'report'
        
    case 'doit'
        % model (convolved design matrix, so after aamod_firstlevel_design)
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        SPM = loadbetter(spmpath);
        % epi data is loaded implicitly by reading off SPM
        nchunks = length(SPM.nscan);

        % get condition names and concatenate runs
        names = arrayfun(@(x)x.name{1},SPM.Sess(1).U,...
            'uniformoutput',false);
        assert(issorted(names),['conditions should be sorted to ' ...
            'avoid scrambling condition order']);
        % take out conditions of no interest
        ts = aap.tasklist.currenttask.settings;
        [validnames,inds] = setdiff(names,ts.ignorelabels);
        ncon = length(validnames);

        % TODO: support for covariates (opt.extraregressors). But I suspect
        % most covariates we can imagine would be subsumed by the noise PCs
        % anyway.

        % convert to GLMdenoise format
        design = cell(1,nchunks);
        epi = cell(1,nchunks);
        frameperiod = SPM.xY.RT;
        interp = [];
        dur = [];
        for chunk = 1:nchunks
            fprintf('converting session %d\n',chunk);
            vols = SPM.xY.P(SPM.Sess(chunk).row,:);
            % TODO: I suspect we will need to use spm_sample_vol to read
            % masked brain here. But we should first assess performance
            % with unmasked.
            epi{chunk} = single(spm_read_vols(spm_vol(vols)));
            % convert onsets in seconds to binary design matrix
            design{chunk} = zeros(SPM.nscan(chunk),ncon);
            % add conditions of interest
            for con = 1:ncon
                conind = inds(con);
                assert(strcmp(names{con},...
                    SPM.Sess(chunk).U(conind).name{1}),...
                    'regressor order must be identical across runs');
                % convert onsets to scans (time0 == volume 1 and so on)
                scantime = 1 + SPM.Sess(chunk).U(conind).ons / frameperiod;
                ons = round(scantime);
                % track interpolation in s
                interp = [interp abs(scantime-ons)*frameperiod];
                assert(all(ons>0 & ons<SPM.nscan(chunk)),...
                    'onsets outside of session scan');
                if isempty(dur)
                    dur = SPM.Sess(chunk).U(conind).dur;
                else
                    assert(all(dur == SPM.Sess(chunk).U(conind).dur),...
                        'all events must have the same duration');
                end
                design{chunk}(ons,con) = 1;
            end
            assert(rank(design{chunk})==size(design{chunk},2),...
                'rank deficient design matrix in chunk %d',chunk);
        end
        fprintf(['interpolated design matrix to scans ' ...
            '(errors [s]: mean=%.2f, min=%.2f, max=%.2f)\n'],...
            mean(interp),min(interp),max(interp));
        clear SPM

        % Kendrick's empirical HRFs look very noisy for short stimuli (<3
        % s) and don't work at all for duration 0 (impulse response). In
        % any case they are extremely similar to the spm_hrf (reassuring!).
        % So we are going to use the SPM HRF when durations are 0. The
        % situation is more complicated for longer durations since you'd
        % then need to convolve the spm_hrf to make a predicted HRF for
        % longer durations.
        % (note that the duration input is basically redundant - Kendrick's
        % code only uses the duration for creating the HRF so when this is
        % given (e.g. when 'optimize' and 'hrfknobs' gives an HRF)
        if dur==0 && isempty(ts.hrfknobs)
            fprintf('event duration 0 so using spm HRF\n');
            ts.hrfknobs = normalizemax(spm_hrf(frameperiod));
        end

        % fit GLMdenoise model
        outdir = fullfile(aas_getsubjpath(aap,subj),'glmdenoise');
        mkdirifneeded(outdir);
        figdir = fullfile(outdir,'diagnostic_figures');
        mkdirifneeded(figdir);
        [results,denoisedepi] = GLMdenoisedata(design,epi,dur,frameperiod,...
            ts.hrfmodel,ts.hrfknobs,ts.opt,figdir);

        % save and describe
        outpath_results = fullfile(outdir,'results.mat');
        save(outpath_results,'results','-v7');
        aap=aas_desc_outputs(aap,subj,'glmdenoise_results',outpath_results);
        outpath_epi = fullfile(outdir,'denoisedepi.mat');
        save(outpath_epi,'denoisedepi','-v7');
        aap=aas_desc_outputs(aap,subj,'glmdenoise_epi',outpath_epi);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
