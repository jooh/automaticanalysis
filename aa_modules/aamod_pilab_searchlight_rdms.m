% generate RDMs for each searchlight in the volume's mask
% [aap,resp]=aamod_pilab_searchlight_rdms(aap,task,subj)
function [aap,resp]=aamod_pilab_searchlight_rdms(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data
        vpath = aas_getfiles_bystream(aap,subj,'pilab_volume');
        vol = loadbetter(vpath);

        % get searchlights
        spherepath = aas_getfiles_bystream(aap,subj,...
            'pilab_searchlight_spheres');
        spheres = loadbetter(spherepath);

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
        data = NaN([npairs vol.nfeatures vol.desc.samples.nunique.chunks]);
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpaths_sessrdms = [];

        % run
        assert(vol.desc.samples.nunique.chunks>0,...
          'vol must have defined chunks in meta.samples');
        for sess = 1:vol.desc.samples.nunique.chunks
            % copying here saves memory per worker in parfor
            sessvol = vol(vol.meta.samples.chunks==sess,:);
            sessdata = data(:,:,sess);
            fprintf('running searchlight %d of %d...\n',sess,...
              vol.desc.samples.nunique.chunks);
            tic;
            parfor n = 1:vol.nfeatures
                % skip empty spheres (these come out as NaN)
                if ~any(spheres.data(n,:))
                    continue
                end
                % pull data direct rather than make instance for speed
                sphdata = sessvol.data(:,spheres.data(n,:));
                sessdata(:,n) = pdist(sphdata,...
                    aap.tasklist.currenttask.settings.distancemetric);
            end
            fprintf('finished in %s.\n',seconds2str(toc));
            % RDMs
            data(:,:,sess) = sessdata;
            % make a volume instance
            sessdisvol = MriVolume(sessdata,vol);
            outpath_sessdata = fullfile(pidir,sprintf(...
                'searchlight_rdms_session%02d.mat',sess));
            save(outpath_sessdata,'sessdisvol');
            outpaths_sessrdms = [outpaths_sessrdms; outpath_sessdata];
        end

        % make average RDM across sessions and save
        meandata = mean(data,3);
        disvol = MriVolume(meandata,vol);
        outpath_mean = fullfile(pidir,'searchlight_rdms_mean.mat');
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
