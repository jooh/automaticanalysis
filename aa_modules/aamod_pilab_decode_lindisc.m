% fit discriminant to data in each ROI (only ROI support at the moment).
%
% [aap,resp]=aamod_pilab_decode_lindisc(aap,task,subj)
function [aap,resp]=aamod_pilab_decode_lindisc(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        designpath = aas_getfiles_bystream(aap,subj,'pilab_design');
        designvol = loadbetter(designpath);
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);
        ts = aap.tasklist.currenttask.settings;

        [~,~,ext] = fileparts(ts.contrastpath);
        if strcmp(lower(ext),'.mat')
            % load
            contrasts = loadbetter(ts.contrastpath);
        else
            % assume script that returns the needed struct
            contrasts = feval(ts.contrastpath);
        end

        % split the data into cell arrays
        [designcell,epicell] = splitvol(ts.split,designvol,epivol);

        % get ROIs / spheres
        roipath = aas_getfiles_bystream(aap,subj,...
            'pilab_rois');
        rois = loadbetter(roipath);

        % check that parfor is available
        if ~matlabpool('size')
            try
                matlabpool;
            catch
                warning('no matlabpool available')
            end
        end

        % prepare output
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        sumdata = [];
        nsplit = length(designcell);
        nanmask = false([nsplit rois.nsamples]);

        if ischar(ts.cvsplit)
            ts.cvsplit = eval(ts.cvsplit);
        end

        % run the beast
        for sp = 1:nsplit
            fprintf('running rois for split %d of %d...\n',sp,nsplit);
            % check for nans
            nanmask = ~any(isnan(epicell{sp}.data),1);
            if ~all(nanmask)
                epicell{sp} = epicell{sp}(:,nanmask);
            end
            % intersect rois and designvol
            allok = epicell{sp}.mask & rois.mask;
            if ~isequal(allok,epicell{sp}.mask)
                epicell{sp} = epicell{sp}(:,epicell{sp}.linind2featind(...
                    find(allok)));
            end
            if ~isequal(allok,rois.mask)
                rois = rois(:,rois.linind2featind(find(allok)));
            end
            [splitres(sp),splitnull(sp),splitboot(sp)] = ...
                roidata_lindisc(rois,designcell{sp},epicell{sp},...
                contrasts,'sgolayK',ts.sgolayK,'sgolayF',ts.sgolayF,...
                'split',ts.cvsplit,'covariatedeg',ts.covariatedeg,...
                'targetlabels',ts.targetlabels,'ignorelabels',...
                ts.ignorelabels,'glmclass',ts.glmclass,'glmvarargs',...
                ts.glmvarargs,'nperm',ts.nperm,'nboot',ts.nboot,...
                'usegpu',ts.usegpu);
            nanmask(sp,:) = any(isnan(splitres(sp).t),1);
        end % sp 1:nsplit

        if ts.usegpu
            splitres = gather(splitres);
            splitnull = gather(splitnull);
            splitboot = gather(splitboot);
            nanmask = gather(nanmask);
        end

        % remove any nan rois from all splits
        anynan = any(nanmask,1);
        if any(anynan)
            for sp = 1:nsplit
                for fi = {'cols_roi','t','medianboot','medianste',...
                        'ppara','pperm','nfeatures'}
                    fstr = fi{1};
                    if isfield(splitres,fstr)
                        splitres(sp).(fstr)(:,anynan,:) = [];
                    end
                    if isfield(splitnull,fstr)
                        splitnull(sp).(fstr)(:,anynan,:) = [];
                    end
                    if isfield(splitboot,fstr)
                        splitboot(sp).(fstr)(:,anynan,:) = [];
                    end
                end
            end
            nnans = sum(anynan);
            fprintf(['removed %d NaN ROIs from analysis ' ...
                '(%.2f%% of total).\n'],nnans,...
                100*(nnans/length(anynan)));
        end

        % make mean result across splits
        meanres = collapsestruct(splitres);
        meannull = collapsestruct(splitnull);
        meanboot = collapsestruct(splitboot);
        % recompute p and boot stats based on distribution of means
        meanres.p_mean = permpvalue(meannull.t);
        [meanres.medianboot_mean,meanres.medianste_mean] = bootprctile(...
            meanboot.t);

        % write out mean results
        % res
        outpath_meandata = fullfile(pidir,'decoder_t_mean.mat');
        save(outpath_meandata,'meanres');
        aap=aas_desc_outputs(aap,subj,'pilab_decoder_t_mean',...
            outpath_meandata);
        % nulldist
        outpath_meannull = fullfile(pidir,'decoder_nulldist_mean.mat');
        save(outpath_meannull,'meannull');
        aap=aas_desc_outputs(aap,subj,'pilab_decoder_nulldist_mean',...
            outpath_meannull);
        % bootdist
        outpath_meanboot = fullfile(pidir,'decoder_bootdist_mean.mat');
        save(outpath_meanboot,'meanboot');
        aap=aas_desc_outputs(aap,subj,'pilab_decoder_bootdist_mean',...
            outpath_meanboot);

        % write out session results
        outpaths_splitres = [];
        outpaths_splitnull = [];
        outpaths_splitboot = [];
        for sp = 1:nsplit
            % res
            outpath_sessdata = fullfile(pidir,sprintf(...
                'decoder_t_split%02d.mat',sp));
            res = splitres(sp);
            save(outpath_sessdata,'res');
            outpaths_splitres = [outpaths_splitres; outpath_sessdata];
            % nulldist
            outpath_nulldata = fullfile(pidir,sprintf(...
                'decoder_nulldist_split%02d.mat',sp));
            nulldist = splitnull(sp);
            save(outpath_nulldata,'nulldist');
            outpaths_splitnull = [outpaths_splitnull; outpath_nulldata];
            % bootdist
            outpath_bootdata = fullfile(pidir,sprintf(...
                'decoder_bootdist_split%02d.mat',sp));
            bootdist = splitboot(sp);
            save(outpath_bootdata,'bootdist');
            outpaths_splitboot = [outpaths_splitboot; outpath_bootdata];
        end
        aap=aas_desc_outputs(aap,subj,'pilab_decoder_t_sess',...
            outpaths_splitres);
        aap=aas_desc_outputs(aap,subj,'pilab_decoder_nulldist_sess',...
            outpaths_splitnull);
        aap=aas_desc_outputs(aap,subj,'pilab_decoder_bootdist_sess',...
            outpaths_splitboot);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
