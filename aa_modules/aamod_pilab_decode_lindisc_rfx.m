% group analysis of discriminant results.
%
% [aap,resp]=aamod_pilab_decode_lindisc_rfx(aap,task)
function [aap,resp]=aamod_pilab_decode_lindisc_rfx(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);
        % save results in main module directory
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        ts = aap.tasklist.currenttask.settings;

        % TODO - outputmode. For now we assume ROI-based. 
        % just load all the results in one go
        for s = 1:nsub
            subres(s) = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_decoder_t_mean'));
        end

        % find all possible ROIs (not all subjects will have all ROIs
        % necessarily)
        allrois = {subres.cols_roi};
        allrois = horzcat(allrois{:});
        urois = unique(allrois);
        nroi = length(urois);

        % but we do assume that contrasts are identical.  we could support
        % this case in theory but let's not for now since it likely
        % indicates a user error.
        concell = {subres.rows_contrast};
        assert(isequal(concell{:}),...
            'different contrasts in different subjects');
        concell = concell{1};
        ncon = length(concell);

        % so now we know what the group result is going to be like
        % (TODO - maybe support other data fields, e.g. median boot)
        dat = NaN([ncon nroi nsub]);
        groupres = struct('rows_contrast',{subres(1).rows_contrast},...
            'cols_roi',{urois},'t',dat,'nfeatures',NaN([1 nroi nsub]));

        % populate the groupres struct
        for s = 1:nsub
            for r = 1:length(subres(s).cols_roi)
                thisroi = subres(s).cols_roi{r};
                indroi = strcmp(groupres.cols_roi,thisroi);
                groupres.t(:,indroi,s) = subres(s).t(:,r);
                groupres.nfeatures(1,indroi,s) = subres(s).nfeatures(1,r);
            end
        end

        % now we can get the mean struct 
        meanres = struct('rows_contrast',{groupres.rows_contrast},...
            'cols_roi',{groupres.cols_roi},'z_subject',...
            {{aap.acq_details.subjects.mriname}});
        meanres.mean = nanmean(groupres.t,3);
        meanres.n = sum(~isnan(groupres.t),3);
        meanres.std = nanstd(groupres.t,[],3);
        meanres.sterr = meanres.std ./ sqrt(meanres.n);
        % nb t is now the group t stat (on the t stats), and ppara / pperm
        % are computed on the subjects too
        meanres.t = meanres.mean ./ meanres.sterr;
        meanres.ppara = tcdf(-meanres.t,meanres.n-1);
        meanres.nfeatures = nanmean(groupres.nfeatures,3);
        meanres.nfeatures_std = nanstd(groupres.nfeatures,[],3);

        % check that sample size is the same for all contrasts by comparing
        % each to the first
        t = bsxfun(@eq,meanres.n,meanres.n(1,:));
        assert(all(t(:)),...
            'mismatched sample size across contrasts for a given ROI');

        % now we want to a) bootstrap, b) permutation test. Hm. The
        % preferred and maximally future-proof version would be to create a
        % GLM with one feature per ROI and run it separately on each
        % contrast (with a constant). That way we can use sample-based
        % permutation and bootstrap methods. The catch is that we have lots
        % of NaNs in the data here so we would probably wind up with a
        % different model for each contrast and ROI. This is a pain but it
        % is the 'correct' way to do things.
        meanres.pperm = NaN(size(meanres.t));
        meanres.medianboot = NaN(size(meanres.t));
        meanres.medianste = NaN(size(meanres.t));

        % for ease of indexing, data is now
        % subject by contrast by roi
        roidata = shiftdim(groupres.t,2);
        if ts.nperm > 1 || ts.nboot > 0
            % we only have to bother with all this fuss if we are doing
            % some kind of resampling
            for r = 1:nroi
                % indexing the first row works because n is the same for
                % each contrast
                design = ones(meanres.n(1,r),1);
                thisroidata = roidata(:,:,r);
                % any and all are assumed to be equivalent due to test
                % above (nb, magically design is already the right length)
                goodroi = ~any(isnan(thisroidata),2);
                model = GLM(design,thisroidata(goodroi,:));

                % permutation test
                if ts.nperm > 1
                    nulldist = permutesamples(model,ts.nperm,'fit',...
                        [1 ncon]);
                    meanres.pperm(:,r) = permpvalue(nulldist);
                end

                % bootstrap
                if ts.nboot > 0
                    bootdist = bootstrapsamples(model,ts.nboot,'fit',...
                        [1 ncon]);
                    [meanres.medianboot(:,r), meanres.medianste(:,r)] = ...
                        bootprctile(bootdist);
                end
            end
        end

        % save and describe
        outpath = fullfile(pidir,'decoder_t_rfx.mat');
        save(outpath,'meanres');
        aap=aas_desc_outputs(aap,'pilab_decoder_t_rfx',...
            outpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
