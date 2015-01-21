% generate RDMs for each ROI (whether it's a set of searchlights or ROIs).
% Generalised version of aamod_pilab_searchlight_rdms.
%
% [aap,resp]=aamod_pilab_rdms(aap,task,subj)
function [aap,resp]=aamod_pilab_rdms(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data
        vpath = aas_getfiles_bystream(aap,subj,'pilab_volume');
        vol = loadbetter(vpath);

        % get ROIs / spheres
        roipath = aas_getfiles_bystream(aap,subj,...
            'pilab_rois');
        rois = loadbetter(roipath);

        % check that parfor is available
        if ~matlabpool('size')
            try
                matlabpool local
            catch
                warning('no matlabpool available')
            end
        end

        % prepare output
        assert(~isempty(vol.desc.samples.nunique.labels),...
          'input vol must have defined labels');
        npairs = nchoosek(vol.desc.samples.nunique.labels,2);
        % npairs by nrois by nchunks
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpaths_sessrdms = [];

        % run
        assert(vol.desc.samples.nunique.chunks>0,...
          'vol must have defined chunks in meta.samples');
        sumdata = zeros([npairs rois.nsamples]);
        % track NaN features - may appear in different runs if nan masking
        nanmask = false([vol.desc.samples.nunique.chunks rois.nsamples]);
        sessdisvolcell = cell(vol.desc.samples.nunique.chunks,1);
        for sess = 1:vol.desc.samples.nunique.chunks
            fprintf('running rois for session %d of %d...\n',sess,...
              vol.desc.samples.nunique.chunks);
            sessdisvolcell{sess} = roidata2rdmvol(rois,vol(...
                vol.meta.samples.chunks==sess,:),...
                aap.tasklist.currenttask.settings.distancemetric);
            nanmask(sess,:) = any(isnan(sessdisvolcell{sess}.data),1);

            sumdata = sumdata + sessdisvolcell{sess}.data;
        end
        % now remove any nan features from all sessdisvols
        anynan = any(nanmask,1);
        if any(anynan)
            nnans = sum(anynan);
            fprintf(['removed %d NaN ROIs from analysis ' ...
                '(%.2f%% of total).\n'],nnans,...
                100*(nnans/length(anynan)));
        end
        sessdisvolcell = cellfun(@(dv)dv(:,~anynan),sessdisvolcell,...
            'uniformoutput',false);
        % and from sums
        sumdata(:,anynan) = [];
        % extract meta features for mean rdm vol (needs to be after main
        % loop to avoid nan ROIs)
        mfeatures = sessdisvolcell{1}.meta.features;
        for sess = 1:vol.desc.samples.nunique.chunks
            if isfield(sessdisvolcell{sess}.meta.features,'nfeatures')
                fn = sprintf('nfeatures_split%02d',sess);
                mfeatures.(fn) = sessdisvolcell{sess}.meta.features.nfeatures;
            end
        end

        % write out sessions
        for sess = 1:vol.desc.samples.nunique.chunks
            outpath_sessdata = fullfile(pidir,sprintf(...
                'rdms_session%02d.mat',sess));
            sessdisvol = sessdisvolcell{sess};
            save(outpath_sessdata,'sessdisvol');
            outpaths_sessrdms = [outpaths_sessrdms; outpath_sessdata];
        end
        % make average RDM across sessions and save
        disvol = SPMVolume(sumdata/vol.desc.samples.nunique.chunks,...
            sessdisvolcell{1},'metafeatures',mfeatures);
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
