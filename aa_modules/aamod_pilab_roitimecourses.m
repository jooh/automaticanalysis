% convert pilab_epi from voxels to ROI time courses for each ROI in
% pilab_rois
%
% [aap,resp]=aamod_pilab_roitimecourses(aap,task,subj)
function [aap,resp]=aamod_pilab_roitimecourses(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);

        % get ROIs / spheres
        roipath = aas_getfiles_bystream(aap,subj,...
            'pilab_rois');
        rois = loadbetter(roipath);

        maskepivol = epi2roitimecourse(rois,epivol);

        save(epipath,'maskepivol')
        aap=aas_desc_outputs(aap,subj,'pilab_epi',epipath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
