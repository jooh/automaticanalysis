% normalise dissimilarities to template brain
% [aap,resp]=aamod_pilab_searchlight_rdms_average(aap,task,subj)
function [aap,resp]=aamod_pilab_rdms_normalise(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get RDMs
        sessrdmpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_sess');
        meanrdmpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
        matname = aas_getfiles_bystream(aap,subj,'normalisation_seg_sn');
        meanrdms = loadbetter(meanrdmpath);
        mpath = aas_getfiles_bystream(aap,'pilab_mask_group');
        % normal, group level mask
        mask = spm_read_vols(spm_vol(mpath)) > 0;

        % mean
        fprintf('normalising mean rdms...')
        meandir = fullfile(fileparts(meanrdmpath),'niftivol_mean');
        mkdirifneeded(meandir);
        meanrdms_norm = dumptodiskandnorm(meanrdms,matname,meandir,mask);
        save(meanrdmpath,'meanrdms_norm');
        fprintf('done.\n');
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_mean',...
            meanrdmpath);

        % sessions
        nsess = size(sessrdmpath,1);
        if nsess==1
            % can just shortcut this since it would produce identical
            % output
            save(sessrdmpath,'meanrdms_norm');
        else
            % otherwise, iterate over sessions...
            for s = 1:nsess
                fprintf(...
                    'normalising rdms for session %d of %d...',s,nsess);
                sessdir = fullfile(fileparts(sessrdmpath(s,:)),sprintf(...
                    'niftivol_session_%03d',s));
                mkdirifneeded(sessdir);
                sessrdms = loadbetter(sessrdmpath(s,:));
                sessrdms_norm = dumptodiskandnorm(sessrdms,matname,...
                    sessdir,mask);
                save(sessrdmpath(s,:),'sessrdms_norm');
                fprintf('done.\n');
            end
        end
        aap=aas_desc_outputs(aap,subj,'pilab_data_rdms_sess',...
            sessrdmpath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

function normdisvol = dumptodiskandnorm(disvol,segpath,outdir,normmask)

% first dump each dissimilarity to disk
filenames = cell(disvol.nsamples,1);
mask = disvol.mask;
V = disvol.header;
datamat = disvol.data;

fprintf('dumping %d dissimilarities to disk...',disvol.nsamples);
tic;
parfor d = 1:disvol.nsamples
    filenames{d} = fullfile(outdir,sprintf('dissimilarity_%04d.nii',d));
    datavec2nifti(datamat(d,:),mask,filenames{d},V);
end

maskind = find(normmask);
nfeat = numel(maskind);
featind = 1:nfeat;
% no preallocation of datavecs because it's parfor
VI = spm_vol(char(filenames));
% normalise in memory (spm_write_sn doesn't return all headers unless you
% wrap it in a loop like this)
fprintf('\nnormalising dissimilarities...');
parfor v = 1:length(VI)
    VO(v) = spm_write_sn(VI(v),segpath);
    datavecs(v,featind) = VO(v).dat(maskind);
end
fprintf('finished in %s\n',seconds2str(toc));

% make new instance with header (minus the data)
normdisvol = SPMVolume(datavecs,normmask,'header',rmfield(VO(1),'dat'));
assert(~any(isnan(normdisvol.data(:))),'nans in normalised disvol');
