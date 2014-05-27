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

        % get stimuli
        spath = aas_getfiles_bystream(aap,subj,'pilab_stimuli');
        stimuli = loadbetter(spath);
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);

        ts = aap.tasklist.currenttask.settings;

        rdms = feval(ts.predictorfun,rdms_org);
        % save and describe
        save(rdmpath,'rdms');

        figdir = fullfile(pidir,'figures');
        mkdirifneeded(figdir);
        aap = aas_desc_outputs(aap,subj,'pilab_rsapredictors',rdmpath);

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
        if all(isnan(rsm(:)))
            fprintf('RSM is all nan, skipping.\n');
        else
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
        end
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
