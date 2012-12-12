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

        % iterate over sessions
        nsess = length(aap.acq_details.selected_sessions);
        % prepare output
        rdms = NaN([vol.nlabels vol.nlabels vol.nfeatures nsess]);
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpaths_sessrdms = [];

        % run
        for sess = 1:nsess
            % copying here saves memory per worker in parfor
            sessvol = vol(vol.chunks==sess,:);
            sessrdms = rdms(:,:,:,sess);
            fprintf('running searchlight %d of %d...\n',sess,nsess)
            tic;
            parfor n = 1:vol.nfeatures
                % skip empty spheres (these come out as NaN)
                if ~any(spheres(n,:))
                    continue
                end
                sphvol = sessvol(:,spheres(n,:));
                sessrdms(:,:,n) = squareform(pdist(sphvol.data,...
                    aap.tasklist.currenttask.settings.distancemetric));
            end
            fprintf('finished in %s.\n',seconds2str(toc));
            % RDMs
            rdms(:,:,:,sess) = sessrdms;
            outpath_rdms = fullfile(pidir,sprintf(...
                'searchlight_rdms_session%02d.mat',sess));
            save(outpath_rdms,'sessrdms');
            outpaths_sessrdms = [outpaths_sessrdms; outpath_rdms];
        end

        % make average RDM across sessions and save
        meanrdms = mean(rdms,4);
        outpath_meanrdm = fullfile(pidir,'searchlight_rdms_mean.mat');
        save(outpath_meanrdm,'meanrdms');

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sess',...
            outpaths_sessrdms);
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',...
            outpath_meanrdm);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
