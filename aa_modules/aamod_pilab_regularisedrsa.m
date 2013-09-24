% Use regularised (ridge) regression to compare a set of model RDMs to
% data RDMs from each split.
%
% [aap,resp]=aamod_pilab_regularisedrsa(aap,task,subj)
function [aap,resp]=aamod_pilab_regularisedrsa(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data RDMs
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_sess');
        ndata = size(vpath,1);
        disvols = arrayfun(@(x)loadbetter(vpath(x,:)),(1:ndata)',...
            'uniformoutput',false);

        % predictor RDMs
        predictpath = aas_getfiles_bystream(aap,subj,...
            'pilab_rsapredictors');
        predictors = loadbetter(predictpath);
        npredictors = length(predictors);

        % check that parfor is available
        if ~matlabpool('size')
            try
                matlabpool local
            catch
                warning('no matlabpool available')
            end
        end

        ts = aap.tasklist.currenttask.settings;

        % parse models
        % rescore 'all' keyword to mean all available predictors
        allinds = find(cellfun(@(it)strcmp(it,'all'),{ts.models.rdminds}));
        for ind = allinds
            ts.models(ind).rdminds = 1:npredictors;
        end
        % add reduced models for each individual predictor
        if ts.eachmodelalone
            for p = 1:npredictors
                ts.models(end+1) = struct('name',predictors(p).name,...
                    'rdminds',p);
            end
        end
        nmodels = length(ts.models);

        % prepare outputs
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);

        % streams
        % result mats with parameter estimates etc
        modelrespaths = {};
        % bootstrap mats
        bootr2paths = {};
        bootestimatepaths = {};
        % niftis of median R2 / estimate - only for searchlights
        r2paths = {};
        estpaths = {};
        % model comparisons now in separate module - so if you just want
        % estimates (e.g. for second level) you can run this with a single
        % model.
        if ischar(ts.ktotry)
            ts.ktotry = eval(ts.ktotry);
        end
        roinames = disvols{1}.meta.features.names;

        % store r2 for each model to compare
        bootr2 = NaN([ts.nboot nmodels disvols{1}.nfeatures]);
        for modind = 1:nmodels
            modname = ts.models(modind).name;
            fprintf('running model %s (%d of %d)...\n',modname,modind,...
                nmodels);
            tic;
            % prepare output for this model
            modresdir = fullfile(resdir,modname);
            mkdirifneeded(modresdir);
            figdir = fullfile(modresdir,'figures');
            mkdirifneeded(figdir);
            result = struct;

            % create model
            glm = disvol2glm(predictors(ts.models(modind).rdminds),...
                disvols,'RidgeRSA');
            dataind = 1:ndata;

            % bootstrap the entire model - crossvalidation and all
            bootmedianr2 = NaN([glm(1).nfeatures,ts.nboot]);
            bootmedianestimates = NaN([glm(1).npredictors ...
                glm(1).nfeatures ts.nboot]);
            % store k parameter for each bootstrap
            bootk = NaN([glm(1).nfeatures ts.nboot]);

            parfor b = 1:ts.nboot
                % draw a random sample of the conditions in the glm with
                % replacement 
                if b == 1
                    % use true model in first boot - also helps getting
                    % estimates when not bootstrapping
                    bootglm = glm;
                else
                    bootglm = glm.bootsampleconditions;
                end
                % leave one out over run splits
                splitr2 = NaN([ndata glm(1).nfeatures]);
                splitestimates = NaN([ndata glm(1).npredictors ...
                    glm(1).nfeatures]);
                splitk = NaN([ndata glm(1).nfeatures]);

                for sp = dataind
                    testind = dataind==sp;
                    trainind = ~testind;
                    % tune regularisation with crossvalidation on the
                    % training runs
                    splitk(sp,:) = bootglm(trainind).crossvalidateproperty(...
                        'k',ts.ktotry);
                    % assign winning k to entire dataset
                    [bootglm.k] = deal(splitk(sp,:));
                    % test prediction performance on the test run
                    splitr2(sp,:) = bootglm(testind).rsquare(...
                        bootglm(trainind).predictY(bootglm(testind).X));
                    % get parameter estimates for test split with chosen k
                    splitestimates(sp,:,:) = bootglm(testind).fit;
                end

                % store median across splits
                bootmedianr2(:,b) = squeeze(median(splitr2,1));
                bootmedianestimates(:,:,b) = squeeze(median(...
                    splitestimates,1));
                bootk(:,b) = squeeze(median(splitk,1));
            end

            % store for model comparison
            bootr2(:,modind,:) = bootmedianr2';

            % get percentile median / st error estimates
            result.mediank = median(bootk,2);
            [result.r2,result.r2err] = bootprctile(bootmedianr2);
            [result.estimates,result.esterr] = bootprctile(...
                bootmedianestimates);
            % save mats of results and bootstrap distributions
            [result.bootdist_r2path,bootr2paths{end+1}] = deal(fullfile(...
                modresdir,'bootdist_r2.mat'));
            save(bootr2paths{end},'bootmedianr2','-v7.3');
            [result.bootdir_estpath,bootestimatepaths{end+1}] = deal(...
                fullfile(modresdir,'bootdist_estimates.mat'));
            save(bootestimatepaths{end},'bootmedianestimates','-v7.3');
            modelrespaths{end+1} = fullfile(modresdir,'results.mat');
            save(modelrespaths{end},'result','-v7.3');
            % diagnostic figures and nifti outputs depend on analysis type
            bnames = {predictors(ts.models(modind).rdminds).name};
            x = 1:glm(1).npredictors;

            switch ts.outputmode
                case 'searchlight'
                    % write out nifti of r2 map
                    r2paths{end+1} = fullfile(modresdir,'r2.nii');
                    disvols{1}.data2file(result.r2,r2paths{end});
                    % and a diagnostic figure
                    F = slicefigure(disvols{1}.data2mat(result.r2),...
                        ts.ylims,'median crossvalidated R^2');
                    printstandard(fullfile(figdir,'r2'));
                    close(F);

                    % write niftis of each individual median estimate
                    for b = x
                        estpaths{end+1} = fullfile(modresdir,sprintf(...
                            'medest%02d-%s.nii',bnames{b}));
                        disvols{1}.data2file(result.estimates(b,:),...
                            estpaths{end});
                    end

                case 'roi'
                    % make a single figure showing r2 for each ROI
                    F = barchart(result.r2,'errors',result.r2err,...
                        'labels',stripbadcharacters(roinames,''));
                    ylabel({'median crossvalidated R^2',...
                        '\pm 1 standard error'});
                    title([modname ' model']);
                    printstandard(fullfile(figdir,'r2'));
                    close(F);
                    if glm(1).npredictors>1
                        % make diagnostic bar chart of estimates + errors
                        % for each ROI
                        for r = 1:disvols{1}.nfeatures
                            roistr = roinames{r};
                            F = barchart(result.estimates(:,r),'errors',...
                                result.esterr(:,r),'labels',...
                                stripbadcharacters(bnames,' '));
                            ylabel({'median parameter estimate',...
                                '\pm 1 standard error'});
                            titlestr = stripbadcharacters(roistr);
                            if isfield(disvols{1}.meta.features,'nfeatures')...
                                    && ~isempty(...
                                    disvols{1}.meta.features.nfeatures)
                                titlestr = sprintf('%s (%d voxels)',...
                                    titlestr,...
                                    disvols{1}.meta.features.nfeatures(r));
                            end
                            title(titlestr);
                            printstandard(fullfile(figdir,sprintf(...
                                'bar_%s',roistr)));
                            close(F);
                        end
                    else
                        % just a single bar chart
                        F = barchart(result.estimates,'errors',...
                            result.esterr,'labels',stripbadcharacters(...
                            roinames,''));
                        title(modname);
                        ylabel({'median parameter estimate',...
                            '\pm 1 standard error'});
                        printstandard(fullfile(figdir,sprintf(...
                            'bar_allrois')));
                        close(F);
                    end
                otherwise
                    error('unknown outputmode: %s',ts.outputmode)
            end
            fprintf('finished %d bootstraps in %s\n',ts.nboot,...
                seconds2str(toc));
        end

        switch ts.outputmode
            case 'roi'
                % plot r2 distributions for each ROI
                for r = 1:disvols{1}.nfeatures
                    roistr = roinames{r};
                    F = figure;
                    P = histline(bootr2(:,:,r),linspace(0,1,40),...
                        'linewidth',2);
                    L = legend(P,{ts.models.name},'location','best');
                    % TODO: add vertical line showing 'true' model
                    % performance for each
                    set(L,'box','off');
                    xlabel('median crossvalidated R^2');
                    ylabel('frequency')
                    titlestr = stripbadcharacters(roistr);
                    if isfield(disvols{1}.meta.features,'nfeatures')...
                            && ~isempty(...
                            disvols{1}.meta.features.nfeatures)
                        titlestr = sprintf('%s (%d voxels)',...
                            titlestr,...
                            disvols{1}.meta.features.nfeatures(r));
                    end
                    title(titlestr);
                    box off
                    printstandard(fullfile(resdir,sprintf(...
                        'modelcomp_%s',roistr)));
                    close(F);
                end
            case 'searchlight'
                % maybe scatter of voxel r2s for each model pair?
                keyboard;
                % TODO
            otherwise
                error('unknown outputmode: %s',ts.outputmode)
        end

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_rsaresult',...
            modelrespaths);
        aap=aas_desc_outputs(aap,subj,'pilab_bootdist_r2',...
            bootr2paths);
        aap=aas_desc_outputs(aap,subj,'pilab_bootdist_estimate',...
            bootestimatepaths);
        if ~isempty(r2paths)
            % niftis - only for searchlight mode
            aap=aas_desc_outputs(aap,subj,'pilab_r2',...
                r2paths);
            aap=aas_desc_outputs(aap,subj,'pilab_estimates',...
                estpaths);
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
