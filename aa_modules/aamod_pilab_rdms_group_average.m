% average RDMs from aamod_pilab_rdms across subjects
%
% [aap,resp]=aamod_pilab_rdms_group_average(aap,task);
function [aap,resp]=aamod_pilab_rdms_group_average(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);
        % save results in main module directory
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        switch aap.tasklist.currenttask.settings.outputmode
            case 'searchlight'
                % relatively simple - we just iterate over subjects and sum
                % up all the rdms, then divide by nsub to get group mean
                for s = 1:nsub
                    meanrdms = loadbetter(aas_getfiles_bystream(aap,s,...
                        'pilab_data_rdms_mean'));
                    sessrdms = cellfun(@(fpath)loadbetter(fpath),...
                        aas_getfiles_bystream(aap,s,...
                        'pilab_data_rdms_sess'),'uniformoutput',false);
                    sesspaths = aas_getfiles_bystream(aap,s,...
                        'pilab_data_rdms_sess');
                    sessrdms = arrayfun(@(x) loadbetter(...
                        sesspaths(x,:)),'uniformoutput',false);
                    if s == 1
                        sumrdms_mean = meanrdms.data;
                        sumrdms_sess = cellfun(@(d)d.data,sessrdms,...
                            'uniformoutput',false);
                        xyz_ref = meanrdms.xyz;
                        nsess = length(sessrdms);
                    else
                        assert(all(xyz_ref(:))==meanrdms.xyz(:),...
                            'mismatched coordinates! need to normalise?');
                        % I guess this isn't necessarily a problem if you
                        % only want to analyse the group mean across
                        % session, but for now let's just crash here to
                        % catch bugs.
                        assert(length(sessrdms)==nsess,...
                            'mismatched session numbers across subjects');
                        sumrdms_mean = sumrdms_mean + meanrdms.data;
                        sumrdms_sess = arrayfun(@(x)sumrdms_sess{x} + ...
                            sessrdms{x}.data,1:length(sessrdms),...
                            'uniformoutput',false);
                    end
                end % s 1:nsub
                grouprdms_mean = sumrdms_mean / nsub;
                grouprdms_sess = cellfun(@(d) d / nsub,sumrdms_sess,...
                    'uniformoutput',false);
                % put in SPMVolume (initialised with meta data from last
                % subject) and save 
                groupdisvol = SPMVolume(grouprdms_mean,meanrdms);
                sessoutpaths = cell(nsess,1);
                for sess = 1:nsess
                    groupdisvol_sess = SPMVolume(grouprdms_sess{sess},...
                        sessrdms{1});
                    sessoutpaths{sess} = fullfile(pidir,sprintf(...
                        'rdms_session_%02d.mat',sess));
                    save(sessoutpaths{sess},'groupdisvol_sess');
                end

            case 'roi'
                % tricky - different subjects may have different ROIs, and
                % we want to preserve meta data on the size of each ROI and
                % how many subjects have it. But no need to use sumrdm
                % approach - probably enough memory to just load all data
                % in one go
                meanrdms_all = {};
                sessrdms_all = {};
                for s = 1:nsub
                    meanrdms_all(end+1,1) = {loadbetter(...
                    aas_getfiles_bystream(aap,s,'pilab_data_rdms_mean'))};
                    sesspaths = aas_getfiles_bystream(aap,s,...
                        'pilab_data_rdms_sess');
                    if s == 1
                        nsess = size(sesspaths,1);
                    else
                        assert(nsess == size(sesspaths,1),...
                            'mismatched session numbers across subjects');
                    end
                    sessrdms_all(end+1,:) = arrayfun(@(x) loadbetter(...
                        sesspaths(x,:)),1:nsess,'uniformoutput',false);
                end
                % get all the possible rois
                allnames = cellfun(@(d)d.meta.features.names(:)',...
                    meanrdms_all,'uniformoutput',false);
                allnames = cat(2,allnames{:});
                % attempt to preserve original ROI order (this will only
                % work if all subjects have all the ROIs)
                unames = unique(allnames,'stable');
                nroi = length(unames);
                fprintf('found %d unique ROIs\n',nroi);
                meanrdms_mat = NaN([meanrdms_all{1}.nsamples nroi nsub]);
                % track number of voxels for meta - I guess it would be
                % ideal to also track other meta data but this will do for
                % now.
                meanrdms_nfeatures = NaN([nroi nsub]);
                % one cell per session
                sessrdms_mat = repmat({meanrdms_mat},[1 nsess]);
                sessrdms_nfeatures = repmat({meanrdms_nfeatures},[1 ...
                nsess]);
                for r = 1:nroi
                    roistr = unames{r};
                    % find subjects with this ROI
                    goodsubs = find(cellfun(@(subd) any(strcmp(...
                        subd.meta.features.names,roistr)),meanrdms_all,...
                        'uniformoutput',true));
                    ngoodsub = numel(goodsubs);
                    fprintf('found %s in %d/%d subjects\n',...
                        roistr,ngoodsub,nsub);
                    for s = 1:ngoodsub
                        goodsind = goodsubs(s);
                        % mean RDM for this subject
                        submean = meanrdms_all{goodsind};
                        hits = strcmp(roistr,submean.meta.features.names);
                        assert(sum(hits)==1,...
                            'expected exactly one matching roi');
                        meanrdms_mat(:,r,goodsind) = submean.data(:,hits);
                        meanrdms_nfeatures(r,goodsind) = ...
                            submean.meta.features.nfeatures(hits);
                        % session RDMs for this subject
                        for sess = 1:nsess
                            subsess = sessrdms_all{goodsind,sess};
                            assert(all(hits == strcmp(roistr,...
                                subsess.meta.features.names)),...
                                'mismatched ROIs in session vs mean data');
                            sessrdms_mat{sess}(:,r,goodsind) = ...
                                subsess.data(:,hits);
                            sessrdms_nfeatures{sess}(r,goodsind) = ...
                                subsess.meta.features.nfeatures(hits);
                        end
                    end % s ngoodsub
                end % r nroi

                % create mean dissimilarities across subjects and save
                meanrdms_mat = nanmean(meanrdms_mat,3);
                % careful - here we can't just plug in meta data from a
                % subject because we may have different features (ROIs) in
                % the aggregate
                groupdisvol = SPMVolume(meanrdms_mat,[],'metafeatures',...
                    struct('names',{unames},'nfeatures',nanmean(...
                    meanrdms_nfeatures,2)','nsub',sum(~isnan(...
                    meanrdms_nfeatures),2)'));
                for sess = 1:nsess
                    groupdisvol_sess = SPMVolume(nanmean(...
                        sessrdms_mat{sess},3),[],'metafeatures',...
                        struct('names',{unames},'nfeatures',nanmean(...
                        sessrdms_nfeatures{sess},2)','nsub',sum(isnan(...
                        sessrdms_nfeatures{sess}),2)'));
                    sessoutpaths{sess} = fullfile(pidir,sprintf(...
                        'rdms_session_%02d.mat',sess));
                    save(sessoutpaths{sess},'groupdisvol_sess');
                end

            otherwise
                error('unknown outputmode: %s',...
                    aap.tasklist.currenttask.settings.outputmode);
        end % switch outputmode

        % save group mean (mean across sessions and subjects)
        outpath = fullfile(pidir,'rdms_groupmean.mat');
        save(outpath,'groupdisvol');
        aap = aas_desc_outputs(aap,'pilab_data_rdms_group_mean',outpath);
        % and describe session-specific output (already saved)
        aap = aas_desc_outputs(aap,'pilab_data_rdms_group_sess',sessoutpaths);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
