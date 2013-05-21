% compare each pilab_data_rdms_mean to each predictor in
% pilab_rsapredictors using permutation testing.
% [aap,resp]=aamod_pilab_rsapermtest(aap,task,subj)
function [aap,resp]=aamod_pilab_rsapermtest(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data RDMs
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_sess');
        nsplit = size(vpath,1);
        % load disvols into a cell array
        disvols = arrayfun(@(x)loadbetter(vpath(x,:)),1:nsplit,...
            'uniformoutput',false);
        roinames = disvols{1}.meta.features.names;
        nroi = disvols{1}.nfeatures;

        [rbyroi,rbysplit] = rdmreproducibility(disvols{:});

        % make outputs
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);
        figdir = fullfile(resdir,'figures');
        mkdirifneeded(figdir);

        ts = aap.tasklist.currenttask.settings;
        switch ts.outputmode
            case 'searchlight'
                % write out nifti
                roiout = fullfile(resdir,'rdmrep.nii');
                disvols{1}.data2file(rbyroi,roiout);
                % make diagnostic figure
                if isempty(ts.ylims)
                    ts.ylims = [0 1];
                end
                F = slicefigure(disvols{1}.data2mat(rbyroi),ts.ylims,...
                    'mean spearman rho');
                title('RDM reproducibility');
                printstandard(fullfile(figdir,'rdmrep'));
                close(F);
            case 'roi'
                % just save mat
                roiout = fullfile(resdir,'rdmrep.mat');
                save(roiout,'rbyroi');
                % make bar graph 
                F = figure;
                x = 1:nroi;
                B = bar(x,rbyroi,.6,'edgecolor','none',...
                    'facecolor',[.6 .6 .6]);
                ylabel('mean spearman rho')
                title('RDM reproducibility');
                set(gca,'xtick',x,'xticklabel',...
                    stripbadcharacters(roinames,' '),'tickdir','out',...
                    'ticklength',get(gca,'ticklength')*.5);
                xlim([x(1)-1 1+x(end)]);
                if ~isempty(ts.ylims)
                    ylim(ts.ylims);
                end
                rotateXLabels(gca,45);
                box off
                printstandard(fullfile(figdir,'rdmrep'));
                close(F);
            otherwise
                error('unrecognised outputmode setting: %s',...
                    ts.outputmode);
        end

        % make split-based plot - same regardless of input data
        splitout = fullfile(resdir,'rdmrep_bysplit.mat');
        save(splitout,'rbysplit');
        F = figure;
        x = 1:nsplit;
        B = bar(x,rbysplit,.6,'edgecolor','none','facecolor',[.6 .6 .6]);
        ylabel('mean spearman rho')
        if ~isempty(ts.ylims)
            ylim(ts.ylims);
        end
        set(gca,'xtick',x);
        xlim([x(1)-1 1+x(end)]);
        xlabel('split');
        box off
        printstandard(fullfile(figdir,'rdmrep_bysplit'));
        close(F);

        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_rdmrep_byroi',roiout);
        aap=aas_desc_outputs(aap,subj,'pilab_rdmrep_bysplit',splitout);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
