% group-level fixed effects inference for RSA.
%
% [aap,resp]=aamod_pilab_rsa_ffx(aap,task)
function [aap,resp]=aamod_pilab_rsa_ffx(aap,task)

resp='';

switch task
    case 'doit'
        ts = aap.tasklist.currenttask.settings;
        nsub = length(aap.acq_details.subjects);

        % get data RDMs
        switch ts.input
            case 'mean'
                disvol = arrayfun(@(x)loadbetter(aas_getfiles_bystream(...
                    aap,x,'pilab_data_rdms_mean')),1:nsub,...
                    'uniformoutput',0);
            case 'sess'
                for s = 1:nsub
                    vpath = aas_getfiles_bystream(aap,s,...
                        'pilab_data_rdms_sess');
                    disvol{s,1} = arrayfun(@(x)loadbetter(vpath(x,:)),...
                        1:size(vpath,1),'uniformoutput',0);
                end
                disvol = cat(1,disvol{:});
            otherwise
                error('unknown input: %s',ts.input);
        end
        % remove ROIs with missing data and ensure everything is in a
        % consistent order (note that unlike the RFX modules less than full
        % sample sizes are not supported here)
        [disvol{1:numel(disvol)}] = intersectvols(disvol{:});

        % predictor RDMs
        if isempty(ts.predictorfun)
            logstr(['using the first subjects'' pilab_rsapredictors ' ...
                'for FFX group analysis\n']);
            predictpath = aas_getfiles_bystream(aap,1,...
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
            'customfun',ts.customfun,...
            'contrasts',ts.contrasts,...
            'customfits',ts.customfits,...
            'bootmeth',ts.bootmeth,...
            'bootprep',ts.bootprep);
        logstr('finished in %s.\n',seconds2str(toc));

        % save and describe
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);

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

        switch ts.outputmode
            case 'roi'
                % save result as mat
                outpath_r = fullfile(pidir,'rsa_r_ffx.mat');
                save(outpath_r,'res');
                aap=aas_desc_outputs(aap,'pilab_res_ffx',...
                    outpath_r);

            case 'searchlight'
                ncon = size(res.r,1);
                rpaths = cell(ncon,1);
                % write out r and p niftis for each predictor
                for c = 1:ncon
                    outname =  stripbadcharacters(res.rows_contrast{c},...
                        '_');
                    rpaths{c} = fullfile(pidir,sprintf(...
                        'rsa_ffx_r_%s.nii',outname));
                    disvol{1}.data2file(res.r(c,:),rpaths{c});
                    % also write out p maps (these are not logged as an
                    % output stream though)
                    if ts.nperm > 1
                        ppath = fullfile(pidir,sprintf(...
                            'rsa_ffx_-log10p_%s.nii',outname));
                        disvol{1}.data2file(-log10(res.pperm(c,:)),ppath);
                        pfwe = permpfwe(squeeze(nulldist.r(c,:,:))');
                        pfwepath = fullfile(pidir,sprintf(...
                            'rsa_ffx_-log10pFWE_%s.nii',outname));
                        disvol{1}.data2file(-log10(pfwe),pfwepath);
                    end
                end
                aap=aas_desc_outputs(aap,'pilab_res_ffx',rpaths);
            otherwise
                error('unknown outputmode: %s',ts.outputmode);
        end

        % null and bootdists get saved as mats regardless
        outpath_null = fullfile(pidir,'rsa_ffx_nulldist.mat');
        save(outpath_null,'nulldist','-v7.3');
        outpath_boot = fullfile(pidir,'rsa_ffx_bootdist.mat');
        save(outpath_boot,'bootdist');

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
