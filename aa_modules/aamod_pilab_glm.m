% generate pilab_result from fitting epivol and
% designvol.
%
% [aap,resp]=aamod_pilab_glm(aap,task,subj)
function [aap,resp]=aamod_pilab_glm(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        designpath = aas_getfiles_bystream(aap,subj,'pilab_design');
        designvol = loadbetter(designpath);
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = loadbetter(epipath);

        ts = aap.tasklist.currenttask.settings;

        args = {'nperm',ts.nperm,'nboot',ts.nboot,...
            'glmclass',ts.glmclass,'glmclassargs',ts.glmclassargs,...
            'customfun',ts.customfun,'customfits',ts.customfits,...
            'contrasts',ts.contrasts,...
            'bootmeth',ts.bootmeth,'bootprep',ts.bootprep,...
            'permmeth',ts.permmeth,'permprep',ts.permprep};
        tic;
        logstr('running %s with %d rois\n',ts.glmclass,epivol.nfeatures);
        [res,nulldist,bootdist] = roidata_glm(designvol,epivol,args{:});
        logstr('finished in %s.\n',seconds2str(toc));

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');

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

        if isempty(ts.outputmode)
            switch class(epivol)
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
                disvol{1}.data2file(res.b(c,:),rpaths);
                % also write out p maps (these are not logged as an
                % output stream though)
                if isfield(res,'pperm')
                    ppath = fullfile(pidir,sprintf(...
                        'glm_-log10p_%s.nii',outname));
                    disvol{1}.data2file(-log10(res.pperm(c,:)),ppath);
                    pfwe = permpfwe(squeeze(nulldist.b(c,:,:))');
                    pfwepath = fullfile(pidir,sprintf(...
                        'glm_-log10pFWE_%s.nii',outname));
                    disvol{1}.data2file(-log10(pfwe),pfwepath);
                end
            end
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
