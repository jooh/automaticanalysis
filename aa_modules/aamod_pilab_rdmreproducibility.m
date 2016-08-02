% get the mean spearman rho for each leave-one-out split of the session
% RDMs.
% [aap,resp]=aamod_pilab_rdmreproducibility(aap,task,subj)
function [aap,resp]=aamod_pilab_rdmreproducibility(aap,task,subj)

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

        % prepare outputs
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);
        figdir = fullfile(resdir,'figures');
        mkdirifneeded(figdir);

        ts = aap.tasklist.currenttask.settings;
        % TODO - support maskfun - re-do analysis with different
        % dissimilarities masked out. Some hairy indexing of disvol rows
        % required.

        ndis = disvols{1}.nsamples;
        ncon = npairs2n(ndis);
        if isempty(ts.maskfun)
            baserdm = ones(ncon,ncon);
            baserdm(logical(eye(ncon))) = 0;
            masks = struct('name','all','RDM',baserdm);
        else
            masks = feval(ts.maskfun);
            % input check
            rdmat = asrdmvec(masks);
            assert(numel(unique(rdmat(~isnan(rdmat))))==1,...
                'rdms must be mean style (one level only)');
        end
        nmask = length(masks);
        for m = 1:nmask
            tic;
            fprintf('running rdmreproducibility analysis %d of %d...',...
                m,nmask);
            % mask out the conditions of no interest
            inds = ~isnan(asrdmvec(masks(m).RDM));
            maskvols = cellfun(@(thisvol)thisvol(inds,:),disvols,...
                'uniformoutput',false);
            [rdmrep{m},rbysplit{m}] = ...
                rdmreproducibility(maskvols{:});
            fprintf(' finished in %s\n',seconds2str(toc));
        end

        switch ts.outputmode
            case 'searchlight'
                for m = 1:nmask
                    % write out nifti
                    roiout{m} = fullfile(resdir,sprintf('rdmrep_%s.nii',...
                        masks(m).name));
                    disvols{1}.data2file(rdmrep{m},roiout{m});
                    % make diagnostic figure
                    if isempty(ts.ylims)
                        ts.ylims = [0 1];
                    end
                    F = slicefigure(...
                        disvols{1}.data2mat(rdmrep{m}),ts.ylims,...
                        'mean spearman rho');
                    title('RDM reproducibility');
                    printstandard(stripextension(roiout{m}));
                    close(F);
                end
            case 'roi'
                error('roi outputmode currently unsupported');
                % just save mat
                roiout = fullfile(resdir,'rdmrep.mat');
                save(roiout,'rdmrep');
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
        for m = 1:nmask
            F = figure;
            x = 1:nsplit;
            B = bar(x,rbysplit{m},.6,'edgecolor','none',...
                'facecolor',[.6 .6 .6]);
            ylabel('mean spearman rho');
            if ~isempty(ts.ylims)
                ylim(ts.ylims);
            end
            set(gca,'xtick',x);
            xlim([x(1)-1 1+x(end)]);
            xlabel('split');
            box off
            printstandard(fullfile(figdir,sprintf('rdmrep_bysplit_%s',...
                masks(m).name)));
            close(F);
            drawnow;
        end

        % describe output - now rsa_r so we can plug into rsa_rfx module
        aap=aas_desc_outputs(aap,subj,'pilab_rsa_r',roiout);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
