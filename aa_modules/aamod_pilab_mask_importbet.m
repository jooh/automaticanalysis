% [aap,resp]=aamod_pilab_importmask_bet(aap,task,subj)
function [aap,resp]=aamod_pilab_importmask_bet(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get mask 
        mpath = aas_getfiles_bystream(aap,subj,'epiBETmask');
        % first mask is the brain mask
        V = spm_vol(mpath(1,:));
        mask = spm_read_vols(V) > 0;

        subdir = aas_getsubjpath(aap,subj);
        pidir = fullfile(subdir,'pilab');
        mkdirifneeded(pidir);
        outpath = fullfile(pidir,'pilab_mask.nii');
        V.fname = outpath;
        spm_write_vol(V,mask);
        aap=aas_desc_outputs(aap,subj,'pilab_mask',outpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
