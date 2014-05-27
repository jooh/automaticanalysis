% Import predictor RDMs for RSA.
%
% [aap,resp]=aamod_pilab_importrsapredictors(aap,task,subj)
function [aap,resp]=aamod_pilab_importrsapredictors(aap,task,subj)

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
        outpath = fullfile(pidir,'pilab_rsapredictors.mat');

        if ~isempty(ts.resortind)
            if ischar(ts.resortind)
                ts.resortind = feval(ts.resortind,numel(stimuli));
            end
            rdmmat = asrdmmat(rdms);
            rdmmat = rdmmat(ts.resortind,ts.resortind,:);
            % assume the stimuli are already sorted appropriately
            for c = 1:length(rdms)
                rdms(c).RDM = rdmmat(:,:,c);
            end
        end
        save(outpath,'rdms');

        figdir = fullfile(pidir,'figures');
        mkdirifneeded(figdir);
        aap = aas_desc_outputs(aap,subj,'pilab_rsapredictors',outpath);

        % Also make quick diagnostic figures
        for r = 1:length(rdms)
            % base rdm
            F = figure;
            rdmplot(gca,rdms(r),'labels',{stimuli.image},'nrows',...
                ts.nrows);
            printstandard(fullfile(figdir,['predictor_rdm_' rdms(r).name]));
            close(F);
            im = intensity2rgb(rdms(r).RDM,cmap_bwr);
            imwrite(im,fullfile(figdir,sprintf('predictor_rawrdm_%s.png',...
                rdms(r).name)),'PNG');
        end
        % show second-order RSM
        rsm = corr(asrdmvec(rdms),'rows','pairwise','type','spearman');
        F = figurebetter([],'huge',1);
        subplot(1,2,1);
        [ax,intmap,cmap] = rdmplot(gca,rsm,'limits',[-1 1]);
        axis(ax,'on');
        set(ax,'yticklabel',stripbadcharacters({rdms.name},' '),...
            'xticklabel',[]);
        sax = subplot(1,2,2);
        c = colorbarbetter(sax,intmap,cmap,'label','spearman rho');
        printstandard(fullfile(figdir,'predictor_so_rsm'));
        close(F);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
