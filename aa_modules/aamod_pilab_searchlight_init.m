% make a 2D nfeat by nfeat logical matrix that provides indices for each
% searchlight sphere.
% [aap,resp]=aamod_pilab_searchlight_init(aap,task,subj)
function [aap,resp]=aamod_pilab_searchlight_init(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data
        % vpath = aas_getfiles_bystream(aap,subj,'pilab_volume');
        % roivol = load(vpath);
        % roivol = roivol.roivol;
        
        % load mask
        mpath = aas_getfiles_bystream(aap,subj,'pilab_mask');
        
        %V = spm_vol(mpath(1,:));
        %mask = spm_read_vols(V) > 0;

        % configure searchlight
        ts = aap.tasklist.currenttask.settings;

        % check that parfor is available
        if ~matlabpool('size')
            warning('no matlabpool available')
        end

        [roivol,diagnostic] = mask2searchrois(mpath(1,:),...
            ts.searchlighttype,ts.searchlightparameter);

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);

        % save spheres
        % now as volume 
        outpath_spheres = fullfile(pidir,'searchlight_spheres.mat');
        save(outpath_spheres,'roivol');

        % save niftis of diagnostics
        outpath_n = fullfile(pidir,...
            'diagnostic_searchlight_nvoxpersphere.nii');
        roivol.data2file(diagnostic.nsphere,outpath_n);
        outpath_s = fullfile(pidir,...
            'diagnostic_searchlight_nspherepervox.nii');
        roivol.data2file(diagnostic.nsampled,outpath_s);
        outpath_r = fullfile(pidir,'diagnostic_searchlight_radius.nii');
        roivol.data2file(diagnostic.r,outpath_r);

        % describe outputs
        % (now treat searchlights like any ROI)
        aap=aas_desc_outputs(aap,subj,'pilab_rois',...
            outpath_spheres);
        aap=aas_desc_outputs(aap,subj,'pilab_searchlight_radius',...
            outpath_r);
        aap=aas_desc_outputs(aap,subj,'pilab_searchlight_nvox',...
            outpath_n);
        aap=aas_desc_outputs(aap,subj,'pilab_searchlight_nspheres',...
            outpath_s);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



