% dummy module to stop processing at a given stage

function [aap,resp]=aamod_crashhere(aap,task,subj)

resp='';

switch task
    case 'doit'
        error('this is where we crash');
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
