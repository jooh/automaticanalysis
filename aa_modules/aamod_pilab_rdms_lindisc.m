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

        % get ROIs / spheres
        roipath = aas_getfiles_bystream(aap,subj,...
            'pilab_rois');
        rois = loadbetter(roipath);

        if ~isempty(ts.setclass)
            fprintf('setting data to %s\n',ts.setclass);
            epivol.data = feval(ts.setclass,epivol.data);
            designvol.data = feval(ts.setclass,designvol.data);
        end
        % make sure we have the same ROIs and voxels across splits
        [rois,epivol] = intersectvols(rois,epivol);
        % now that the ROIs and voxels are in register this should reduce
        % memory use considerably
        validvox = any(rois.data~=0,1);
        % at the moment we assume this for detecting searchlight maps. So
        % need to crash if this is not met
        if rois.nsamples == rois.nfeatures
            assert(all(validvox),'unused voxels in searchlight map');
        else
            rois = rois(:,validvox);
            epivol = epivol(:,validvox);
        end

        % check that parfor is available
        if ~matlabpool('size')
            warning('no matlabpool available')
        end

        % prepare output
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');

        [disvol,nulldist,splitdisvolcell,sessnulldist] = ...
            roidata2rdmvol_lindisc_batch(rois,...
            designvol,epivol,'split',ts.split,'glmvarargs',...
            ts.glmvarargs,'cvsplit',ts.cvsplit,'glmclass',ts.glmclass,...
            'sterrunits',ts.sterrunits,'crossvalidate',ts.crossvalidate,...
            'crosscon',ts.crosscon,'searchvol',ts.searchvol,...
            'demean',ts.demean,'nperms',ts.nperms);

        outpath_mean = fullfile(pidir,'rdms_mean.mat');
        save(outpath_mean,'disvol');

        % save session disvols
        outpaths_sessrdms = [];
        for sp = 1:length(splitdisvolcell)
            outpath_sessdata = fullfile(pidir,sprintf(...
                'rdms_split%02d.mat',sp));
            sessdisvol = splitdisvolcell{sp};
            save(outpath_sessdata,'sessdisvol');
            outpaths_sessrdms = [outpaths_sessrdms; outpath_sessdata];
        end

        outpath_null = fullfile(pidir,'rdms_nulldist.mat');
        save(outpath_null,'nulldist');
        outpath_sessnull = fullfile(pidir,'rdms_sessnulldist.mat');
        save(outpath_sessnull,'sessnulldist');

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_nulldist',...
            outpath_null);
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sessnulldist',...
            outpath_sessnull);
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sess',...
            outpaths_sessrdms);
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',...
            outpath_mean);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
