% compare each pilab_data_rdms_mean to each predictor in
% pilab_rsapredictors using permutation testing.
% [aap,resp]=aamod_pilab_rsapermtest(aap,task,subj)
function [aap,resp]=aamod_pilab_rsapermtest(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data RDMs
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
        disvol = loadbetter(vpath);

        % predictor RDMs
        predictpath = aas_getfiles_bystream(aap,subj,...
            'pilab_rsapredictors');
        predictors = loadbetter(predictpath);

        ts = aap.tasklist.currenttask.settings;
        [rpaths,ppaths,pfwepaths,nulldistpaths] = rsapermtest_batch(...
            disvol,predictors,'outputmode',ts.outputmode,...
            'nperms',ts.nperms,'resdir',fullfile(aas_getsubjpath(aap,...
            subj),'pilab','results'));

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_r',...
            rpaths);
        aap=aas_desc_outputs(aap,subj,'pilab_p',...
            ppaths);
        aap=aas_desc_outputs(aap,subj,'pilab_p_fwe',...
            pfwepaths);
        aap=aas_desc_outputs(aap,subj,'pilab_nulldist',...
            nulldistpaths);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
