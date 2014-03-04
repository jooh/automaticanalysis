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


        % parse input arguments
        if ~iscell(ts.glmvarargs)
            if isempty(ts.glmvarargs)
                ts.glmvarargs = {};
            else
                ts.glmvarargs = {ts.glmvarargs};
            end
        end
        % there are 3 split levels here: 
        % supersplit: defines which part of the data is used for
        % searchlight mapping and which is used to generate the final
        % independent RDM
        % split: defines which sessions are used for computing RDMs (the
        % final RDM for the searchlight / ROI is the average over all
        % splits)
        % cvsplit: defines how the LDt RDM is cross-validated within each
        % session
        if ischar(ts.cvsplit)
            ts.cvsplit = eval(ts.cvsplit);
        end
        if ischar(ts.supersplit)
            ts.supersplit = eval(ts.supersplit);
        end
        if isempty(ts.testsplit)
            % same split in train and test (only works if the split
            % contains the same number of runs in train and test)
            ts.testsplit = ts.split;
        end
        usuper = unique(ts.supersplit);
        nsuper = length(usuper);
        uchunks = epivol.desc.samples.unique.chunks;
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
        predictors(ts.removepredictorinds) = [];
        ncon = length(predictors);
        nsizes = length(ts.ns);

        % for each super-split (roi defining and testing split)
        for super = 1:nsuper
            fprintf('running super split %d of %d...\n',super,nsuper)
            % figure out which data will be used to identify ROI and
            % which to test
            thissuper = usuper(super);

            % find roi masks (possibly split-specific)
            roidir = fullfile(ts.roiroot,...
                aap.acq_details.subjects(subj).mriname);
            if ~isempty(ts.subdir)
                if iscell(ts.subdir)
                    roidir = fullfile(roidir,ts.subdir{super});
                else
                    roidir = fullfile(roidir,ts.subdir);
                end
            end
            assert(exist(roidir,'dir')~=0,'no roi dir found: %s',roidir);
            fprintf('loading rois from directory %s\n',roidir);
            maskvol = roidir2vol(roidir);
            % make sure we have the same masks, ROIs and voxels across splits
            [rois,epivol,maskvol] = intersectvols(rois,epivol,maskvol);

            % for each mask
            for mask = 1:maskvol.nsamples
                maskname = maskvol.meta.samples.names{mask};
                fprintf('\nrunning %s (%d of %d)...\n',maskname,mask,...
                    maskvol.nsamples);
                maskind = maskvol.data(mask,:)~=0;
                epimask = epivol(:,maskind);
                roimask = rois(maskind,maskind);

                fprintf('1. defining searchlight rois...\n')
                % compute searchlight RDMs in training split (separately
                % for different sessions if desired)
                trainind = ts.supersplit~=thissuper;
                trainchunks = uchunks(trainind);
                epitrain = epimask.selectbymeta('chunks',trainchunks);
                designtrain = designvol.selectbymeta('chunks',trainchunks);
                traindisvol = roidata2rdmvol_lindisc_batch(roimask,...
                    designtrain,epitrain,'split',ts.split,'glmvarargs',...
                    ts.glmvarargs,'cvsplit',ts.cvsplit,'glmclass',...
                    ts.glmclass,'sterrunits',ts.sterrunits,...
                    'crossvalidate',ts.crossvalidate,...
                    'batchsize',ts.batchsize);
                % get some RSA effect estimate for each searchlight RDM in
                % the training split
                trainres = roidata_rsa(traindisvol,predictors,...
                    'nperm',1,'nboot',0,'rsaclass',ts.rsaclass,...
                    'rsaclassargs',ts.rsaclassargs,...
                    'splitrsaclass',ts.splitrsaclass,...
                    'splitrsaclassargs',ts.splitrsaclassargs);

                % make a new roivol containing ROIs generating according to
                % each con / ns combination
                roimat = false([ncon*nsizes epimask.nfeatures]);
                definingcon = cell(ncon*nsizes,1);
                names = definingcon;
                thresholds = NaN([ncon*nsizes,1]);
                nvec = thresholds;
                nr = 0;
                for con = 1:ncon
                    for n = 1:nsizes
                        nr = nr+1;
                        thissize = ts.ns(n);
                        % make the ROI volume in the test data according to
                        % the ranked effect in train data
                        [roimat(nr,:),thresholds(nr)] = ...
                            selectbysearchlight(roimask,trainres.r(con,:),...
                            thissize,'selectmode','union');
                        nvec(nr) = thissize;
                        definingcon{nr} = trainres.rows_contrast{con};
                        names{nr} = sprintf('%s %s (%d voxels)',...
                            maskname,definingcon{nr},thissize);
                    end % n = 1:nsizes
                end % con = 1:ncon

                roispheres.(maskname){super} = MriVolume(roimat,epimask,...
                    'metasamples',struct('names',{names},...
                    'thresholds',thresholds,'definingcon',{definingcon},...
                    'nfeatures_target',nvec,'masknames',...
                    {repmat({maskname},[nr 1])}));
                % done with training data - clear here to be extra safe
                % with independence
                clear epitrain designtrain trainind trainchunks traindisvol

                % TEST DATA PROCESSING
                fprintf('2. calculating RDMs for new rois\n');
                testind = ts.supersplit==thissuper;
                testchunks = uchunks(testind);
                epitest = epimask.selectbymeta('chunks',testchunks);
                designtest = designvol.selectbymeta('chunks',testchunks);
                % and make the RDMs (we ignore the session RDMs here for
                % simplicity)
                splitdisvol{mask,super} = roidata2rdmvol_lindisc_batch(...
                    roispheres.(maskname){super},designtest,epitest,...
                    'split',ts.testsplit,'glmvarargs',ts.glmvarargs,...
                    'cvsplit',ts.cvsplit,'glmclass',ts.glmclass,...
                    'sterrunits',ts.sterrunits,...
                    'crossvalidate',ts.crossvalidate,...
                    'batchsize',ts.batchsize);
                % clear for independence safety
                clear epitest designtest testind testchunks
            end % mask = 1:maskvol.nsamples 

            fprintf('---------------------------------\n');
        end % super = 1:nsuper

        % combine the splits into one disvol
        for mask = 1:maskvol.nsamples
            for super = 1:nsuper
                disdata{super} = splitdisvol{mask,super}.data;
                disfeat(super) = splitdisvol{mask,super}.meta.features;
            end
            meandis = mean(cat(3,disdata{:}),3);
            % combine the meta data from all splits for roi and disvol
            roimeta = collapsestruct(disfeat,@mean);
            maskdisvol{mask} = BaseVolume(meandis,...
                'metafeatures',roimeta);
        end
        % collapse all the mask ROIs to one big volume
        % (here we may run into mask trouble)
        meandisvol = cat(2,maskdisvol{:});
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpath_mean = fullfile(pidir,'rdms_mean.mat');
        save(outpath_mean,'meandisvol');
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',...
            outpath_mean);

        % save session disvols
        outpaths_sessrdms = [];
        for super = 1:nsuper
            sessdisvol = cat(2,splitdisvol{:,super});
            outpath_sessdata = fullfile(pidir,sprintf(...
                'rdms_split%02d.mat',super));
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
