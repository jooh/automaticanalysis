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
            contrasts = structdata2class(contrasts,class(epivol.data));
        end

        % get ROIs / spheres
        roipath = aas_getfiles_bystream(aap,subj,...
            'pilab_rois');
        rois = loadbetter(roipath);

        % make sure we have the same ROIs and voxels across splits
        [rois,epivol] = intersectvols(rois,epivol);
        % now that the ROIs and voxels are in register this should reduce
        % memory use considerably
        validvox = any(rois.data~=0,1);
        rois = rois(:,validvox);
        epivol = epivol(:,validvox);

        % split the data into cell arrays
        [designcell,epicell] = splitvol(ts.split,designvol,epivol);

        % prepare output
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        sumdata = [];
        nsplit = length(designcell);

        if ischar(ts.cvsplit)
            ts.cvsplit = eval(ts.cvsplit);
        end

        % run the beast
        for sp = 1:nsplit
            fprintf('running rois for split %d of %d...\n',sp,nsplit);
            tic;
            [splitres(sp),splitnull(sp),splitboot(sp)] = ...
                roidata_lindisc(rois,designcell{sp},epicell{sp},...
                contrasts,'split',ts.cvsplit,'glmclass',ts.glmclass,...
                'glmvarargs',ts.glmvarargs,'nperm',ts.nperm,...
                'nboot',ts.nboot);
            fprintf('finished in %s.\n',seconds2str(toc));
        end % sp 1:nsplit

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
