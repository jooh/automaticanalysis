% generate RDMs for each searchlight in the volume's mask
% [aap,resp]=aamod_pilab_searchlight_rdms(aap,task,subj)
function [aap,resp]=aamod_pilab_searchlight_rdms(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data
        vpath = aas_getfiles_bystream(aap,subj,'pilab_volume');
        vol = load(vpath);
        vol = vol.vol;

        % get searchlights
        spherepath = aas_getfiles_bystream(aap,subj,...
            'pilab_searchlight_spheres');
        spheres = load(spherepath);
        spheres = spheres.spheres;

        % check that parfor is available
        if ~matlabpool('size')
            try
                matlabpool local
            catch
                warning('no matlabpool available')
            end
        end

        % prepare output
        npairs = nchoosek(vol.nlabels,2);
        data = NaN([npairs vol.nfeatures vol.nchunks]);
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpaths_sessrdms = [];

        % run
        for sess = 1:vol.nchunks
            % copying here saves memory per worker in parfor
            sessvol = vol(vol.chunks==sess,:);
            sessdata = data(:,:,sess);
            fprintf('running searchlight %d of %d...\n',sess,vol.nchunks)
            tic;
            parfor n = 1:vol.nfeatures
                % skip empty spheres (these come out as NaN)
                if ~any(spheres(n,:))
                    continue
                end
                sphvol = sessvol(:,spheres(n,:));
                sessdata(:,n) = pdist(sphvol.data,...
                    aap.tasklist.currenttask.settings.distancemetric);
            end
            fprintf('finished in %s.\n',seconds2str(toc));
            % RDMs
            data(:,:,sess) = sessdata;
            % make a volume instance
            sessdisvol = Volume(sessdata,vol);
            outpath_sessdata = fullfile(pidir,sprintf(...
                'searchlight_rdms_session%02d.mat',sess));
            save(outpath_sessdata,'sessdisvol');
            outpaths_sessrdms = [outpaths_sessrdms; outpath_sessdata];
        end

        % make average RDM across sessions and save
        meandata = mean(data,3);
        disvol = Volume(meandata,vol);
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
