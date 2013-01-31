% Restrict analysis to voxels where we obtain a sensible searchlight size
% (default 10) for searchlight diagnostic output. Note that this
% module makes little sense if you have done nvox-based mapping (although I
% guess an equivalent function for that case would be to skip searchlights
% where the radius ends up enormous).
% [aap,resp]=aamod_pilab_restrictbysearchlightn(aap,task,subj)
function [aap,resp]=aamod_pilab_restrictbysearchlightn(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        % mask
        mpath = aas_getfiles_bystream(aap,subj,'freesurfer_gmmask');
        % second is EPI
        V = spm_vol(mpath(2,:));
        mask = spm_read_vols(V);
        % get volume
        vpath = aas_getfiles_bystream(aap,subj,'pilab_volume');
        vol = loadbetter(vpath);
        % searchlight diagnostic
        spath = aas_getfiles_bystream(aap,subj,'pilab_searchlight_nvox');
        xyz = spm_read_vols(spm_vol(spath));
        % intersect to generate new mask
        mask = (mask>0) & (xyz>=aap.tasklist.currenttask.settings.minvox);
        ngone = vol.nfeatures-sum(mask(:)>0);
        fprintf('eliminated %d features (%.2f%% of total)\n',...
          ngone,100*(ngone/vol.nfeatures));
        spm_write_vol(V,mask);
        aap = aas_desc_outputs(aap,subj,'freesurfer_gmmask',mpath);
        % update the volume
        goodind = vol.linind2featind(find(mask));
        vol = vol(:,goodind);
        save(vpath,'vol')
        aap = aas_desc_outputs(aap,subj,'pilab_volume',vpath);
        % and spheres...
        spath = aas_getfiles_bystream(aap,subj,...
            'pilab_searchlight_spheres');
        spheres = loadbetter(spath);
        spheres = spheres(goodind,goodind);
        save(spath,'spheres');
        aap = aas_desc_outputs(aap,subj,'pilab_searchlight_spheres',spath);
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
