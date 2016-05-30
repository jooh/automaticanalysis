% Convert EPI to noise
%
% [aap,resp]=aamod_pilab_epi2noise(aap,task,subj)
function [aap,resp]=aamod_pilab_epi2noise(aap,task,subj)

resp='';

switch task
    case 'doit'
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);
        ts = aap.tasklist.currenttask.settings;
        m = mean(epivol.data);
        sd = std(epivol.data);
        dataorg = epivol.data;
        % convert to gaussian noise with the same mean and standard
        % deviation as the original data
        filterbychunk(epivol,@(x)bsxfun(@times,...
            bsxfun(@plus,randn(size(x)),mean(x)),std(x)));

        % output
        logstr('saving...\n');
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        out = fullfile(pidir,'epivol.mat');
        save(out,'epivol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_epi',out);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
