% AA module - estimate firstlevel_spm model
% [aap,resp]=aamod_firstlevel_estimate(aap,task,subj)
function [aap,resp]=aamod_firstlevel_estimate(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        %get subject directory
        cwd=pwd;
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        % if we want to cut scans outside a timewindow of first/last event
        % for independence purposes
        if ~isinf(aap.tasklist.currenttask.settings.padevents)
            keyboard;
        end
        anadir = fileparts(spmpath);
        cd(anadir)
        SPMdes = spm_fmri_spm_ui(SPM);
        % now check real covariates and nuisance variables are
        % specified correctly
        %% (doesn't seem to be used in any other code)
        %SPMdes.xX.iG=cols_nuisance;
        %SPMdes.xX.iC=cols_interest;
        
        % (maybe disable auto mask, set explicit mask?) Apparently must be
        % done after spm_fmri_spm_ui
        spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
        SPMest = spm_spm(SPMdes);
        cd(cwd);
        % Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));
        
        %  firstlevel_betas (includes related statistical files)
        allbetas=dir(fullfile(anadir,'beta_*'));
        betafns=[];
        for betaind=1:length(allbetas);
            betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
        end
        otherfiles={'mask.hdr','mask.img','ResMS.hdr','ResMS.img','RPV.hdr','RPV.img'};
        for otherind=1:length(otherfiles)
            betafns=strvcat(betafns,fullfile(anadir,otherfiles{otherind}));
        end
        aap=aas_desc_outputs(aap,subj,'firstlevel_betas',betafns);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
