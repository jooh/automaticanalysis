% Create predictor RDMs for idanddistinct study.
% [aap,resp]=aamod_pilab_rsapredictors_oneback(aap,task,subj)
function [aap,resp]=aamod_pilab_rsapredictors_oneback(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;

        % get stimuli
        spath = aas_getfiles_bystream(aap,subj,'pilab_stimuli');
        stimuli = loadbetter(spath);

        ts = aap.tasklist.currenttask.settings;

        rdms = feval(ts.predictorfun,subname,stimuli);

        % save and describe
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);
        outpath = fullfile(pidir,'idanddistinct_rsapredictors.mat');
        save(outpath,'rdms');

        figdir = fullfile(pidir,'figures');
        mkdirifneeded(figdir);
        aap = aas_desc_outputs(aap,subj,'pilab_rsapredictors',outpath);

        % Also make quick diagnostic figures
        for r = 1:length(rdms)
            % base rdm
            F = plotrdms(rdms(r),'labels',{stimuli.image},'nrows',...
                ts.nrows,'titles',rdms(r).name);
            printstandard(fullfile(figdir,['predictor_rdm_' rdms(r).name]));
            close(F);
            im = intensity2rgb(rdms(r).RDM,jet(1e3));
            imwrite(im,fullfile(figdir,sprintf('predictor_rawrdm_%s.png',...
                rdms(r).name)),'PNG');
        end
        % show second-order RSM
        F = figure;
        % as big as it will go...
        set(F,'units','pixels','position',get(0,'screensize'));
        rsm = corr(asrdmvec(rdms),'rows','pairwise','type','spearman');
        F = plotrdms(rsm,'labels',stripbadcharacters({rdms.name},' '),...
            'docb',true,'colorbarargs',{'label','spearman rho'},...
            'titles','second-order similarity','fighand',F,...
            'rotatelabels',90,'cmap',jet(1e3));
        printstandard(fullfile(figdir,'predictor_so_rsm'));
        close(F);
        % mds - nah, too many nans
        %F = plotrdms(1-rsm,'labels',{rdms.name},'domds',true,'titles',...
            %'second-order similarity','dordm',false);
        %printstandard(fullfile(figdir,'predictor_so_mds'));
        %close(F);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
