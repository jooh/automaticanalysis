% Make 2 pilab volume instances: one for the GLM specification, and another
% for its corresponding EPI volumes

function [aap,resp]=aamod_pilab_epiandmodel2volume(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='construct EPI and design pilab Volumes';
        
    case 'summary'
        
    case 'report'
        
    case 'doit'
        mpath = aas_getfiles_bystream(aap,subj,'epiBETmask');
        % first mask is the brain mask
        V = spm_vol(mpath(1,:));
        orgmask = spm_read_vols(V) ~= 0;
        
        % model
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');

        [epivol,designvol] = spm2vol(spmpath,'mask',orgmask,...
            'assumeconvolved=1');

        ts = aap.tasklist.currenttask.settings;
        if aap.tasklist.currenttask.settings.maskbybright
            % replicate GLMdenoise EPI intensity mask
            meanepi = mean(epivol.data,1);
            thresh = prctile(meanepi,ts.brainthresh(1))*ts.brainthresh(2);
            bright = meanepi > thresh;
            % update volume
            epivol = epivol(:,bright);
        end
        % and now the correct mask is
        mask = epivol.mask;
        norg = sum(orgmask(:));
        nnow = sum(mask(:));
        fprintf(['removed %d zero/nan/dark features from epivol ' ...
            'and mask (%.2f%% of total)\n'],norg-nnow,100*...
            ((norg-nnow)/norg));

        % save and describe
        outdir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(outdir);
        % epi
        outpath_epi = fullfile(outdir,'epivol.mat');
        save(outpath_epi,'epivol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_epi',outpath_epi);
        % model
        outpath_design = fullfile(outdir,'designvol.mat');
        save(outpath_design,'designvol','-v7');
        aap=aas_desc_outputs(aap,subj,'pilab_design',outpath_design);
        % the updated mask
        spm_write_vol(V,mask);
        aap=aas_desc_outputs(aap,subj,'epiBETmask',mpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
