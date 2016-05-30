% Preprocess an epivol / designvol pair.
%
% [aap,resp]=aamod_pilab_preprocess(aap,task,subj)
function [aap,resp]=aamod_pilab_preprocess(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        designpath = aas_getfiles_bystream(aap,subj,'pilab_design');
        designvol = loadbetter(designpath);
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);
        ts = aap.tasklist.currenttask.settings;

        [epivol,designvol] = preprocessvols(epivol,designvol,...
            'matchn',ts.matchn,'covariatedeg',ts.covariatedeg,...
            'domedianfilter',ts.medianfilter,'sgdetrend',ts.sgdetrend,...
            'sgolayK',ts.sgolayK,'sgolayF',ts.sgolayF,...
            'dozscore',ts.zscore,'targetlabels',ts.targetlabels,...
            'ignorelabels',ts.ignorelabels,'setclass',ts.setclass,...
            'resortind',ts.resortind,'percentsignal',ts.percentsignal);

        % output
        logstr('saving...\n');
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        out = fullfile(pidir,'epivol.mat');
        save(out,'epivol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_epi',out);
        outdesign = fullfile(pidir,'designvol.mat');
        save(outdesign,'designvol','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_design',outdesign);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
