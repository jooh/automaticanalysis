% compare each pilab_data_rdms_mean to each predictor in
% pilab_rsapredictors.
% [aap,resp]=aamod_pilab_rsa(aap,task,subj)
function [aap,resp]=aamod_pilab_rsa(aap,task,subj)

resp='';

switch task
    case 'doit'
        ts = aap.tasklist.currenttask.settings;

        % get data RDMs
        switch ts.input
            case 'mean'
                vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
                disvol = {loadbetter(vpath)};
            case 'sess'
                vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_sess');
                disvol = arrayfun(@(x)loadbetter(vpath(x,:)),...
                    1:size(vpath,1),'uniformoutput',0);
            otherwise
                error('unknown input: %s',ts.input);
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
        
        % partial rho support. The beauty of this approach is that
        % roidata_rsa does not need to know that this is the analysis we're
        % running...
        if ~isempty(ts.partialpredictors)
            [~,partialind] = intersect({predictors.name},...
                ts.partialpredictors);
            assert(~isempty(partialind),'no match for partialpredictors');
            assert(isempty(ts.rsaclassargs),...
                'rsaclassargs must be empty for partial RSA');
            assert(isempty(ts.splitrsaclass),...
                'no split RSA support for partial RSA at present');
            ts.rsaclassargs = {predictors(partialind)};
            % doesn't make sense to fit these anymore
            predictors(partialind) = [];
        end
        % multirsa and related supports
        if ~isempty(ts.fitmethodpredictors)
            [~,fmind] = intersect({predictors.name},...
                ts.fitmethodpredictors);
            assert(~isempty(fmind),'no match for fitmethodpredictors');
            assert(isempty(ts.fitmethodargs),...
                'fitmethodargs must be empty for this mode');
            assert(isempty(ts.splitrsaclass),...
                'no split RSA support for fitmethodargs at present');
            ts.fitmethodargs = {predictors(fmind)};
            % unlike partialpredictors we still keep the original predictor
            logstr('using %s as fitmethodpredictors\n',...
                ts.fitmethodpredictors);
        end
        ncon = length(predictors);

        logstr('running roidata_rsa with %d predictors and %d rois\n',...
            ncon,disvol{1}.nfeatures);
        tic;
        [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,'nperm',...
            ts.nperm,'nboot',ts.nboot,'rsaclass',ts.rsaclass,...
            'rsaclassargs',ts.rsaclassargs,...
            'splitrsaclass',ts.splitrsaclass,...
            'splitrsaclassargs',ts.splitrsaclassargs,...
            'fitmethod',ts.fitmethod,...
            'fitmethodargs',ts.fitmethodargs,...
            'customfun',ts.customfun);
        logstr('finished in %s.\n',seconds2str(toc));

        % save and describe
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');

        if isempty(ts.outputmode)
            switch class(disvol{1})
                case 'BaseVolume'
                    ts.outputmode = 'roi';
                    logstr('auto-set outputmode to roi\n');
                case 'MriVolume'
                    ts.outputmode = 'searchlight';
                    logstr('auto-set outputmode to searchlight\n');
                otherwise
                    error('unknown outputmode for class %s',class(disvol{1}));
            end
        end

        switch ts.outputmode
            case 'roi'
                % save result as mat
                outpath_r = fullfile(pidir,'rsa_r.mat');
                save(outpath_r,'res');
                aap=aas_desc_outputs(aap,subj,'pilab_rsa_r',...
                    outpath_r);

            case 'searchlight'
                ncon = size(res.r,1);
                rpaths = cell(ncon,1);
                % write out r and p niftis for each predictor
                for c = 1:ncon
                    outname =  stripbadcharacters(res.rows_contrast{c},...
                        '_');
                    rpaths{c} = fullfile(pidir,sprintf('rsa_r_%s.nii',...
                        outname));
                    disvol{1}.data2file(res.r(c,:),rpaths{c});
                    % also write out p maps (these are not logged as an
                    % output stream though)
                    if ts.nperm > 1
                        ppath = fullfile(pidir,sprintf(...
                            'rsa_-log10p_%s.nii',outname));
                        disvol{1}.data2file(-log10(res.pperm(c,:)),ppath);
                        pfwe = permpfwe(squeeze(res.nulldist(c,:,:))');
                        pfwepath = fullfile(pidir,sprintf(...
                            'rsa_-log10pFWE_%s.nii',outname));
                        disvol{1}.data2file(-log10(pfwe),pfwepath);
                    end
                end
                aap=aas_desc_outputs(aap,subj,'pilab_rsa_r',rpaths);
            otherwise
                error('unknown outputmode: %s',ts.outputmode);
        end

        % null and bootdists get saved as mats regardless
        outpath_null = fullfile(pidir,'rsa_nulldist.mat');
        save(outpath_null,'nulldist');
        aap=aas_desc_outputs(aap,subj,'pilab_rsa_nulldist',...
            outpath_null);

        outpath_boot = fullfile(pidir,'rsa_bootdist.mat');
        save(outpath_boot,'bootdist');
        aap=aas_desc_outputs(aap,subj,'pilab_rsa_bootdist',...
            outpath_boot);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
