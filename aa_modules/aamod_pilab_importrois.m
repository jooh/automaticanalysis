% Import roitools ROIs pilab analysis 
function [aap,resp]=aamod_pilab_importrois(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;
        roidir = fullfile(aap.tasklist.currenttask.settings.roiroot,...
            subname);
        assert(exist(roidir,'dir')~=0,'no roi dir found: %s',roidir);

        if aap.tasklist.currenttask.settings.usesmoothrois
            target = 'roivol.mat';
        else
            target = 'roivol_unsm.mat';
        end
        roivol = loadbetter(fullfile(roidir,target));

        % intersect ROI with pilab mask
        %pipath = aas_getfiles_bystream(aap,subj,'pilab_volume');
        %pivol = loadbetter(pipath);
        %goodinds = find(roivol.mask & pivol.mask);
        %roivol = roivol(:,goodinds);

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);
        outpath = fullfile(pidir,'pilab_rois.mat');
        save(outpath,'roivol');
        aap = aas_desc_outputs(aap,subj,'pilab_rois',outpath);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
