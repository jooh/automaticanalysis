% [aap,resp]= aamod_pilab_rdms_searchselect(aap,task,subj)
function [aap,resp]= aamod_pilab_rdms_searchselect(aap,task,subj)

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

        % get ROIs / spheres
        roipath = aas_getfiles_bystream(aap,subj,...
            'pilab_rois');
        rois = loadbetter(roipath);

        % sort out the training predictors
        if isempty(ts.predictorfun)
            predictpath = aas_getfiles_bystream(aap,subj,...
                'pilab_rsapredictors');
            predictors = loadbetter(predictpath);
        else
            fprintf('using custom predictor RDMs from %s\n',...
                ts.predictorfun);
            predictors = feval(ts.predictorfun);
        end
        if ~isempty(ts.selectpredictorinds)
            assert(isempty(ts.removepredictorinds),...
                'cannot both select and remove predictor inds');
            predictors = predictors(ts.selectpredictorinds);
        end
        predictors(ts.removepredictorinds) = [];
        % however you did that, there should now be only a single model RDM
        assert(numel(predictors)==1,'1 predictor only');

        % sort out masks
        % find roi masks (possibly split-specific)
        roidir = fullfile(ts.roiroot,...
            aap.acq_details.subjects(subj).mriname);
        if isempty(ts.subdir)
            % simples
            maskvol = roidir2vol(roidir);
        else
            if ~iscell(ts.subdir)
                maskvol = roidir2vol(fullfile(roidir,ts.subdir));
            else
                for s = 1:numel(ts.subdir)
                    roidir = fullfile(roidir,ts.subdir{super});
                    maskvol{s} = roidir2vol(roidir);
                end
            end
        end

        [meandisvol,sessdisvolcell,roispheres] = ...
            roidata2rdmvol_lindisc_searchselect(rois,designvol,epivol,...
            maskvol,...
            'predictor',predictors,'split',ts.split,'cvsplit',ts.cvsplit,...
            'subsplit',ts.subsplit,'minvoxeln',ts.minvoxeln,...
            'glmclass',ts.glmclass,'glmvarargs',ts.glmvarargs,...
            'sterrunits',ts.sterrunits,...
            'crossvalidate',ts.crossvalidate,...
            'crosscon',ts.crosscon,'maskns',...
            ts.maskns,'rsaclass',ts.rsaclass,'rsaclassargs',...
            ts.rsaclassargs);

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpath_mean = fullfile(pidir,'rdms_mean.mat');
        save(outpath_mean,'meandisvol');
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',...
            outpath_mean);

        % save session disvols
        outpaths_sessrdms = [];
        for sp = 1:length(sessdisvolcell)
            sessdisvol = sessdisvolcell{sp};
            outpath_sessdata = fullfile(pidir,sprintf(...
                'rdms_split%02d.mat',sp));
            save(outpath_sessdata,'sessdisvol');
            outpaths_sessrdms = [outpaths_sessrdms; outpath_sessdata];
        end
        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sess',...
            outpaths_sessrdms);

        % also save ROIs for diagnostic purposes
        sphereout = fullfile(pidir,'searchspheres.mat');
        save(sphereout,'roispheres');
        aap=aas_desc_outputs(aap,subj,'pilab_rois_searchselect',...
            sphereout);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
