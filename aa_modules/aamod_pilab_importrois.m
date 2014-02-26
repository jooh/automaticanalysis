% Import roitools ROIs pilab analysis 
function [aap,resp]=aamod_pilab_importrois(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;
        ts = aap.tasklist.currenttask.settings;
        roidir = fullfile(ts.roiroot,subname);
        if ~isempty(ts.subdir)
            roidir = fullfile(roidir,ts.subdir);
        end
        assert(exist(roidir,'dir')~=0,'no roi dir found: %s',roidir);

        if aap.tasklist.currenttask.settings.usesmoothrois
            target = 'roivol.mat';
        else
            target = 'roivol_unsm.mat';
        end
        potvol = fullfile(roidir,target);
        if exist(potvol,'file') ~= 0
            fprintf('loading rois from roivol %s\n',potvol);
            roivol = loadbetter(fullfile(roidir,target));
        else
            fprintf('loading rois from directory %s\n',roidir);
            roivol = roidir2vol(roidir);
        end

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);
        outpath = fullfile(pidir,'pilab_rois.mat');
        save(outpath,'roivol');
        aap = aas_desc_outputs(aap,subj,'pilab_rois',outpath);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
