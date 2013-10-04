% compare each pilab_data_rdms_mean_group to each predictor in
% pilab_rsapredictors using permutation testing.
% [aap,resp]=aamod_pilab_rsapermtest_group(aap,task)
function [aap,resp]=aamod_pilab_rsapermtest_group(aap,task)

resp='';

switch task
    case 'doit'
        % get data RDMs - group mean here
        vpath = aas_getfiles_bystream(aap,'pilab_data_rdms_group_mean');
        disvol = loadbetter(vpath);

        % predictor RDMs
        % NB! we use the first subject's predictors under the assumption
        % that all subjects have the same predictors. If this is not the
        % case we need to use RFX instead.
        predictpath = aas_getfiles_bystream(aap,1,'pilab_rsapredictors');
        predictors = loadbetter(predictpath);
        fprintf(...
            'loaded predictors for group analysis from first subject\n');

        ts = aap.tasklist.currenttask.settings;
        [rpaths,ppaths,pfwepaths,nulldistpaths] = rsapermtest_batch(...
            disvol,predictors,'outputmode',ts.outputmode,...
            'nperms',ts.nperms,'resdir',fullfile(aas_getstudypath(aap),...
            'pilab','results'));

        % describe outputs
        aap=aas_desc_outputs(aap,'pilab_r_group',...
            rpaths);
        aap=aas_desc_outputs(aap,'pilab_p_group',...
            ppaths);
        aap=aas_desc_outputs(aap,'pilab_p_fwe_group',...
            pfwepaths);
        aap=aas_desc_outputs(aap,'pilab_nulldist_group',...
            nulldistpaths);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
