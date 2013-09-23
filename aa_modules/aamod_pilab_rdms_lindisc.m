% generate RDMs for each ROI (whether it's a set of searchlights or ROIs).
% this variant is different from the standard rdms in that it uses
% pilab_design and pilab_epi directly to fit linear discriminants.
%
% [aap,resp]=aamod_pilab_rdms_lindisc(aap,task,subj)
function [aap,resp]=aamod_pilab_rdms_lindisc(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        designpath = aas_getfiles_bystream(aap,subj,'pilab_design');
        designvol = loadbetter(designpath);
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);
        ts = aap.tasklist.currenttask.settings;

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
        % track NaN features across splits may appear in different runs if
        % nan masking - nans should only really happen here if you enter
        % empty ROIs
        nsplit = length(designcell);
        nanmask = false([nsplit rois.nsamples]);
        splitdiscvolcell = cell(nsplit,1);

        % run the beast
        for sp = 1:nsplit
            fprintf('running rois for split %d of %d...\n',sp,nsplit);
            % cart off to new function
            splitdisvolcell{sp} = roidata2rdmvol_lindisc(rois,...
                designcell{sp},epicell{sp},...
                'sgolayK',ts.sgolayK,'sgolayF',ts.sgolayF,'split',...
                ts.cvsplit,'covariatedeg',ts.covariatedeg,...
                'targetlabels',ts.targetlabels,...
                'ignorelabels',ts.ignorelabels,'glmclass',ts.glmclass,...
                'glmvarargs',ts.glmvarargs,'sterrunits',ts.sterrunits,...
                'crossvalidate',ts.crossvalidate);
            nanmask(sp,:) = any(isnan(splitdisvolcell{sp}.data),1);
            if isempty(sumdata)
                sumdata = splitdisvolcell{sp}.data;
            else
                sumdata = sumdata + splitdisvolcell{sp}.data;
            end
        end % sp 1:nsplit

        % remove any nan features from all splits
        anynan = any(nanmask,1);
        if any(anynan)
            nnans = sum(anynan);
            fprintf(['removed %d NaN ROIs from analysis ' ...
                '(%.2f%% of total).\n'],nnans,...
                100*(nnans/length(anynan)));
        end
        splitdisvolcell = cellfun(@(dv)dv(:,~anynan),splitdisvolcell,...
            'uniformoutput',false);
        % and from sums
        sumdata(:,anynan) = [];

        % extract meta features for mean rdm vol (needs to be after main
        % loop to avoid nan ROIs) and write out
        outpaths_sessrdms = [];
        mfeatures = splitdisvolcell{1}.meta.features;
        for sp = 1:nsplit
            if isfield(splitdisvolcell{sp}.meta.features,'nfeatures')
                fn = sprintf('nfeatures_split%02d',sp);
                mfeatures.(fn) = ...
                    splitdisvolcell{sp}.meta.features.nfeatures;
            end
            outpath_sessdata = fullfile(pidir,sprintf(...
                'rdms_split%02d.mat',sp));
            sessdisvol = splitdisvolcell{sp};
            save(outpath_sessdata,'sessdisvol');
            outpaths_sessrdms = [outpaths_sessrdms; outpath_sessdata];
        end

        % make average RDM across sessions and save
        disvol = MriVolume(sumdata/nsplit,splitdisvolcell{1},...
            'metafeatures',mfeatures);
        outpath_mean = fullfile(pidir,'rdms_mean.mat');
        save(outpath_mean,'disvol');

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sess',...
            outpaths_sessrdms);
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',...
            outpath_mean);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
