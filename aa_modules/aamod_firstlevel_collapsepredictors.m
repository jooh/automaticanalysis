% AA module - collapse predictors in SPM.mat
% [aap,resp]=aamod_firstlevel_collapsepredictors(aap,task,subj)
function [aap,resp]=aamod_firstlevel_collapsepredictors(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        %get subject directory
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        ts = aap.tasklist.currenttask.settings;
        SPM = loadbetter(spmpath);
        SPM = feval(ts.collapsefun,SPM);
        % Describe outputs
        %  firstlevel_spm
        save(spmpath,'SPM');
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',spmpath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
