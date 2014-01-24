% compare each pilab_data_rdms_mean to each predictor in
% pilab_rsapredictors.
% [aap,resp]=aamod_pilab_rsa(aap,task,subj)
function [aap,resp]=aamod_pilab_rsa(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data RDMs
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
        disvol = loadbetter(vpath);

        % predictor RDMs
        ts = aap.tasklist.currenttask.settings;
        predictpath = aas_getfiles_bystream(aap,subj,...
            'pilab_rsapredictors');
        predictors = loadbetter(predictpath);
        predictors(ts.removepredictorinds) = [];
        ncon = length(predictors);


        % this is a bit hacky but sometimes it is advantageous to drop the
        % precision before carting off to roidata_rsa
        if ~isempty(ts.setclass)
            fprintf('setting data and predictors to %s\n',ts.setclass);
            disvol.data = feval(ts.setclass,disvol.data);
            for p = 1:ncon
                predictors(p).RDM = feval(ts.setclass,predictors(p).RDM);
            end
        end

        % cart off to refactored rsapermtest function
        fprintf('running roidata_rsa with %d predictors and %d rois\n',...
            ncon,disvol.nfeatures);
        tic;
        [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,'nperm',...
            ts.nperm,'nboot',ts.nboot,'rsaclass',ts.rsaclass,...
            'rsaclassargs',ts.rsaclassargs,...
            'splitrsaclass',ts.splitrsaclass,...
            'splitrsaclassargs',ts.splitrsaclassargs);
        fprintf('finished in %s.\n',seconds2str(toc));

        % save and describe
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');

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
                    disvol.data2file(res.r(c,:),rpaths{c});
                    % also write out p maps (these are not logged as an
                    % output stream though)
                    if ts.nperm > 1
                        ppath = fullfile(pidir,sprintf(...
                            'rsa_-log10p_%s.nii',outname));
                        disvol.data2file(-log10(res.pperm(c,:)),ppath);
                        pfwe = maxstatpfwe(squeeze(res.nulldist(c,:,:))');
                        pfwepath = fullfile(pidir,sprintf(...
                            'rsa_-log10pFWE_%s.nii',outname));
                        disvol.data2file(-log10(pfwe),pfwepath);
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
