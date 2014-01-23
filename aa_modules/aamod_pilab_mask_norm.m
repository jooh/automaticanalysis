% [aap,resp]=aamod_pilab_mask_norm(aap,task)
function [aap,resp]=aamod_pilab_mask_norm(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);

        for s = 1:nsub
            % get mask 
            mpath = aas_getfiles_bystream(aap,s,'pilab_mask');
            matname = aas_getfiles_bystream(aap,s,'normalisation_seg_sn');
            % normalise 
            VO(s) = spm_write_sn(mpath,matname);
        end

        % binarise
        % this mask will be slightly expanded due to smoothing
        xyz = spm_read_vols(VO) > 0;
        % and intersect across subjects - only voxels that remain in all
        % subjects  after normalisation are counted
        groupmask = all(xyz,4);

        % write out the group mask
        VO = VO(1);
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        VO.fname = fullfile(pidir,'pilab_mask_group.nii');
        spm_write_vol(VO,groupmask);

        % describe 
        aap=aas_desc_outputs(aap,'pilab_mask_group',VO.fname);

        % make a diagnostic figure
        F = figure;
        np = ceil(sqrt(nsub+1));
        for s = 1:nsub
            subplot(np,np,s);
            % 2 where subject overlaps groupmask, 1 where subject mask is
            % getting cropped (subject but no group). The converse case
            % (group but not subject) is impossible since we intersect
            % rather than take the union.
            imagesc(makeimagestack(double(xyz(:,:,:,s)) + ...
                double(groupmask)));
            axis equal;
            axis off;
            title(sprintf('%d: %s',s,aap.acq_details.subjects(s).mriname));
        end
        subplot(np,np,nsub+1);
        imagesc(makeimagestack(groupmask));
        axis equal;
        axis off;
        title('group intersection');
        printstandard(fullfile(pidir,'diagnostic_mask_group'));
        close(F);
        drawnow;

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
