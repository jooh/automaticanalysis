% Restrict analysis to voxels where we obtain a sensible searchlight size
% defined in nvox and/or radius.
% [aap,resp]=aamod_pilab_restrictbysearchlightsize(aap,task,subj)
function [aap,resp]=aamod_pilab_restrictbysearchlightsize(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        % searchlight diagnostic
        spath = aas_getfiles_bystream(aap,subj,'pilab_searchlight_nvox');
        xyz_n = spm_read_vols(spm_vol(spath));
        rpath = aas_getfiles_bystream(aap,subj,'pilab_searchlight_radius');
        xyz_r = spm_read_vols(spm_vol(rpath));
                % and spheres...
        spath = aas_getfiles_bystream(aap,subj,...
            'pilab_rois');
        spheres = loadbetter(spath);
        % intersect to generate new mask
        ts = aap.tasklist.currenttask.settings;
        mask = spheres.mask;
        % pilab mask
        mpath = aas_getfiles_bystream(aap,subj,'pilab_searchlight_nvox');
        mV = spm_vol(mpath);
        mxyz = spm_read_vols(mV) ~= 0;
        assert(isequal(mxyz,mask),...
            'mismatched pilab_mask and searchlight mask')
        mask = (mask>0) & (xyz_r >= ts.minradius) & ...
            (xyz_r <= ts.maxradius) & (xyz_n >= ts.minvox) & ...
            (xyz_n <= ts.maxvox);
        ngone = spheres.nfeatures-sum(mask(:));
        fprintf('eliminated %d features (%.2f%% of total)\n',...
          ngone,100*(ngone/spheres.nfeatures));
        mind = find(mask);
        goodind = ismember(spheres.linind,mind);
        spheres = spheres(goodind,goodind);
        save(spath,'spheres');
        aap = aas_desc_outputs(aap,subj,'pilab_rois',spath);
        spm_write_vol(mV,mask);
        aap = aas_desc_outputs(aap,subj,'pilab_mask',mpath);
        % aaand diagnostics
        for dia = {'nvox','radius','nspheres'}
            streamname = ['pilab_searchlight_' dia{1}];
            dpath = aas_getfiles_bystream(aap,subj,streamname);
            dV = spm_vol(dpath);
            dxyz = spm_read_vols(dV);
            dxyz(~mask) = 0;
            spm_write_vol(dV,dxyz);
            aas_desc_outputs(aap,subj,streamname,dpath);
        end
end
