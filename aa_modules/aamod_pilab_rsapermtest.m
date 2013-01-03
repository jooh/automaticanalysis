% compare each pilab_data_rdms_mean to each predictor in
% pilab_rsapredictors using permutation testing.
% [aap,resp]=aamod_pilab_rsapermtest(aap,task,subj)
function [aap,resp]=aamod_pilab_rsapermtest(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get data RDMs
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
        disvol = load(vpath);
        disvol = disvol.disvol;
        % TODO eventually sensitive to multiple vols in cell array here for ROI
        % RSA

        % predictor RDMs
        predictpath = aas_getfiles_bystream(aap,subj,...
            'pilab_rsapredictors');
        predictors = load(predictpath);
        predictors = predictors.rdms;
        npredictors = length(predictors);

        % check that parfor is available
        if ~matlabpool('size')
            try
                matlabpool local
            catch
                warning('no matlabpool available')
            end
        end

        % make outputs
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);
        rpaths = cell(npredictors,1);
        ppaths = cell(npredictors,1);
        pfwepaths = cell(npredictors,1);
        nulldistpaths = cell(npredictors,1);
        % run separate for each predictor to conserve memory (the nulldists
        % end up too big to save in memory otherwise)
        for pre = 1:npredictors
            fprintf('testing %s (%d of %d)...\n',predictors(pre).name,...
                pre,npredictors);
            % permutation test
            tic;
            [r,p,nulldists] = rsapermtest(predictors(pre),...
                disvol.data,aap.tasklist.currenttask.settings.nperms);
            fprintf('finished in %s.\n',seconds2str(toc));
            % obtain Nichols / Holmes-style FWE-corrected p values
            pfwe = maxstatpfwe(nulldists);
            % r map
            rout = fullfile(resdir,sprintf('%s_r.nii',...
                predictors(pre).name));
            disvol.data2file(r,rout);
            rpaths{pre} = rout;
            % log10 p map
            pout = fullfile(resdir,sprintf('%s_-log10p.nii',...
                predictors(pre).name));
            disvol.data2file(-log10(p),pout);
            ppaths{pre} = pout;
            % FWE-corrected p map
            pfweout = fullfile(resdir,sprintf('%s_-log10pFWE.nii',...
                predictors(pre).name));
            disvol.data2file(-log10(pfwe),pfweout);
            pfwepaths{pre} = pfweout;
            % null distributions
            nullout = fullfile(resdir,sprintf('%s_nulldist.mat',...
                predictors(pre).name));
            % save as volume with massive ndata
            nullvol = Volume(nulldists,disvol);
            % for mysterious reasons Matlab cannot save this in any older
            % version
            save(nullout,'nullvol','-v7.3');
            nulldistpaths{pre} = nullout;
            % not much faith in Matlab garbage collection
            clear nullvol
        end
        % describe outputs
        aap=aas_desc_outputs(aap,subj,'pilab_r',...
            rpaths);
        aap=aas_desc_outputs(aap,subj,'pilab_p',...
            ppaths);
        aap=aas_desc_outputs(aap,subj,'pilab_p_fwe',...
            pfwepaths);
        aap=aas_desc_outputs(aap,subj,'pilab_nulldist',...
            nulldistpaths);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
