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
                npath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_nulldist');
                nullrdm = loadbetter(npath);
            case 'sess'
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

        logstr('running %s/%s with %d predictors and %d rois\n',...
            ts.rsaclass, ts.splitrsaclass, ncon,disvol{1}.nfeatures);
        args = {'nperm',ts.nperm,'nboot',ts.nboot,...
            'rsaclass',ts.rsaclass,'rsaclassargs',ts.rsaclassargs,...
            'splitrsaclass',ts.splitrsaclass,...
            'splitrsaclassargs',ts.splitrsaclassargs,...
            'fitmethod',ts.fitmethod,'fitmethodargs',ts.fitmethodargs,...
            'customfun',ts.customfun,...
            'bootmeth',ts.bootmeth,'bootprep',ts.bootprep,...
            'permdis',nullrdm,'customfits',ts.customfits,...
            'contrasts',ts.contrasts};
        tic;
        [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,args{:});
        logstr('finished in %s.\n',seconds2str(toc));

        % save and describe
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');

        if isempty(ts.outputmode)
            switch class(disvol{1})
                case 'Volume'
                    ts.outputmode = 'roi';
                    logstr('auto-set outputmode to roi\n');
                case 'SPMVolume'
                    ts.outputmode = 'searchlight';
                    logstr('auto-set outputmode to searchlight\n');
                otherwise
                    error('unknown outputmode for class %s',class(disvol{1}));
            end
        end

        % save results as mat regardless
        outpath_r = fullfile(pidir,'rsa_r.mat');
        save(outpath_r,'res','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_result',...
            outpath_r);
        outpath_null = fullfile(pidir,'rsa_nulldist.mat');
        save(outpath_null,'nulldist','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_result_nulldist',...
            outpath_null);
        outpath_boot = fullfile(pidir,'rsa_bootdist.mat');
        save(outpath_boot,'bootdist','-v7.3');
        aap=aas_desc_outputs(aap,subj,'pilab_result_bootdist',...
            outpath_boot);

        % maybe also write out niftis
        switch ts.outputmode
            case 'roi'
                % do nothing
            case 'searchlight'
                ncon = size(res.r,1);
                rpaths = cell(ncon,1);
                % write out r and p niftis for each predictor
                for c = 1:ncon
                    outname =  stripbadcharacters(res.rows_contrast{c},...
                        '_');
                    rpaths = fullfile(pidir,sprintf('rsa_r_%s.nii',...
                        outname));
                    disvol{1}.data2file(res.r(c,:),rpaths);
                    % also write out p maps (these are not logged as an
                    % output stream though)
                    if isfield(res,'pperm')
                        ppath = fullfile(pidir,sprintf(...
                            'rsa_-log10p_%s.nii',outname));
                        disvol{1}.data2file(-log10(res.pperm(c,:)),ppath);
                        pfwe = permpfwe(squeeze(nulldist.r(c,:,:))');
                        pfwepath = fullfile(pidir,sprintf(...
                            'rsa_-log10pFWE_%s.nii',outname));
                        disvol{1}.data2file(-log10(pfwe),pfwepath);
                    end
                end
            otherwise
                error('unknown outputmode: %s',ts.outputmode);
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
