% [aap,resp]=aamod_pilab_selectrsapredictors(aap,task,subj)
function [aap,resp]=aamod_pilab_selectrsapredictors(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;

        rdmpath = aas_getfiles_bystream(aap,subj,'pilab_rsapredictors');
        rdms_org = loadbetter(rdmpath);

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);

        ts = aap.tasklist.currenttask.settings;

        rdms = feval(ts.predictorfun,rdms_org);
        % save and describe
        save(rdmpath,'rdms');
        aap = aas_desc_outputs(aap,subj,'pilab_rsapredictors',rdmpath);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
