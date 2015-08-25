% generate pilab_result from RSA GLM.
%
% [aap,resp]=aamod_pilab_rsaglm(aap,task,subj)
function [aap,resp]=aamod_pilab_rsaglm(aap,task,subj)

resp='';

switch task
    case 'doit'
        ts = aap.tasklist.currenttask.settings;

        % get data RDMs
        switch ts.input
            case 'mean'
                vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
                disvol = loadbetter(vpath);
                npath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_nulldist');
                nullrdm = loadbetter(npath);
            case 'sess'
                error('sess input not supported at present');
                vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_sess');
                disvol = arrayfun(@(x)loadbetter(vpath(x,:)),...
                    1:size(vpath,1),'uniformoutput',0);
                npath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_sessnulldist');
                nullrdm = arrayfun(@(x)loadbetter(npath(x,:)),...
                    1:size(npath,1),'uniformoutput',0);
            otherwise
                error('unknown input: %s',ts.input);
        end

        if ts.nperm > 1 && strcmp(ts.permtype,'rsa')
            % remove nullrdm to prevent roidata_rsa from running
            % RDM-based permutation test
            nullrdm = [];
        end

        % predictor RDMs
        if isempty(ts.predictorfun)
            predictpath = aas_getfiles_bystream(aap,subj,...
                'pilab_rsapredictors');
            predictors = loadbetter(predictpath);
        else
            logstr('using custom predictor RDMs from %s\n',...
                ts.predictorfun);
            predictors = feval(ts.predictorfun);
        end
        if ~isempty(ts.selectpredictorinds)
            assert(isempty(ts.removepredictorinds),...
                'cannot both select and remove predictor inds');
            if ischar(ts.selectpredictorinds)
                ts.selectpredictorinds = eval(ts.selectpredictorinds);
            end
            predictors = predictors(ts.selectpredictorinds);
        end
        predictors(ts.removepredictorinds) = [];
        ncon = length(predictors);

        args = {'nperm',ts.nperm,'nboot',ts.nboot,...
            'glmclass',ts.rsaclass,'glmclassargs',ts.rsaclassargs,...
            'customfun',ts.customfun,'customfits',ts.customfits,...
            'contrasts',ts.contrasts,...
            'bootmeth',ts.bootmeth,'bootprep',ts.bootprep,...
            'permmeth',ts.permmeth,'permprep',ts.permprep};
        tic;
        logstr('running %s with %d rois\n',ts.rsaclass,...
            disvol.nfeatures);
        [res,nulldist,bootdist] = roidata_glm(predictors,disvol,args{:});
        logstr('finished in %s.\n',seconds2str(toc));

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');

        if isempty(ts.outputmode)
            switch class(disvol)
                case 'Volume'
                    ts.outputmode = 'roi';
                    logstr('auto-set outputmode to roi\n');
                case 'SPMVolume'
                    ts.outputmode = 'searchlight';
                    logstr('auto-set outputmode to searchlight\n');
                otherwise
                    error('unknown outputmode for class %s',class(disvol));
            end
        end

        % save results as mat regardless
        outpath_r = fullfile(pidir,'glm_b.mat');
        save(outpath_r,'res','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_result',...
            outpath_r);
        outpath_null = fullfile(pidir,'glm_nulldist.mat');
        save(outpath_null,'nulldist','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_result_nulldist',...
            outpath_null);
        outpath_boot = fullfile(pidir,'glm_bootdist.mat');
        save(outpath_boot,'bootdist','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_result_bootdist',...
            outpath_boot);

        % maybe also write out niftis
        if strcmp(ts.outputmode,'searchlight')
            ncon = size(res.b,1);
            rpaths = cell(ncon,1);
            % write out b and p niftis for each predictor
            for c = 1:ncon
                outname =  stripbadcharacters(res.rows_contrast{c},...
                    '_');
                rpaths = fullfile(pidir,sprintf('glm_b_%s.nii',...
                    outname));
                disvol.data2file(res.b(c,:),rpaths);
                % also write out p maps (these are not logged as an
                % output stream though)
                if isfield(res,'pperm')
                    ppath = fullfile(pidir,sprintf(...
                        'glm_-log10p_%s.nii',outname));
                    disvol.data2file(-log10(res.pperm(c,:)),ppath);
                    pfwe = permpfwe(squeeze(nulldist.b(c,:,:)));
                    pfwepath = fullfile(pidir,sprintf(...
                        'glm_-log10pFWE_%s.nii',outname));
                    disvol.data2file(-log10(pfwe),pfwepath);
                end
            end
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
