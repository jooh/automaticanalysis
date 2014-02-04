% group analysis of searchlight rsa results.
%
% [aap,resp]=aamod_pilab_rsa_rfx_searchlight(aap,task)
function [aap,resp]=aamod_pilab_rsa_rfx_searchlight(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);
        % save results in main module directory
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        ts = aap.tasklist.currenttask.settings;

        % group-level mask
        V = spm_vol(aas_getfiles_bystream(aap,'pilab_mask_group'));
        mask = spm_read_vols(V) > 0;

        % need to convert searchlight niftis to subres format
        for s = 1:nsub
            niftis = aas_getfiles_bystream(aap,s,'pilab_rsa_r');
            subres(s) = nifti2roidata(niftis,mask,'r');
        end

        names = {aap.acq_details.subjects.mriname};
        [subres.name] = names{:};

        fprintf('running roidata_rfx with %d subjects \n',nsub);
        tic;
        [meanres,nulldist,bootdist] = roidata_rfx(subres,'nperm',...
            ts.nperm,'nboot',ts.nboot,...
            'targetfield','r','transfun',ts.transfun,'assumeregister',true);
        fprintf('finished in %s.\n',seconds2str(toc));

        % save mats
        save(fullfile(pidir,'meanres.mat'),'meanres');
        save(fullfile(pidir,'nulldist.mat'),'nulldist','-v7.3');
        save(fullfile(pidir,'bootdist.mat'),'bootdist','-v7.3');

        % TODO consider implementing variance smoothing - certainly looks
        % like the mean searchlight map is crazy smooth but the variance
        % less so.

        % save and describe
        % need to dump out each r and p as a nifti
        % probably easiest to first make an MriVolume instance and then use
        % its methods
        resvol = MriVolume(meanres.mean,mask,'header',V);
        ncon = length(meanres.rows_contrast);
        rpaths = cell(ncon,1);
        for con = 1:length(meanres.rows_contrast)
            % descriptive: mean/stdev
            rpaths{con} = fullfile(pidir,sprintf('rfx_mean_%s.nii',...
                meanres.rows_contrast{con}));
            resvol.data2file(resvol.data(con,:),rpaths{con});
            resvol.data2file(meanres.std(con,:),fullfile(pidir,...
                sprintf('rfx_std_%s.nii',meanres.rows_contrast{con})));
            % parametric p value
            ppara = fullfile(pidir,sprintf(...
                'rfx_-log10ppara_%s.nii',...
                meanres.rows_contrast{con}));
            resvol.data2file(-log10(meanres.ppara(con,:)),ppara);
            % bonferroni p
            pparafwe = fullfile(pidir,sprintf(...
                'rfx_-log10pFWEpara_%s.nii',...
                meanres.rows_contrast{con}));
            resvol.data2file(-log10(meanres.ppara(con,:)*sum(mask(:))),...
                pparafwe);

            % maybe perm p values too?
            if ts.nperm > 1
                pperm = fullfile(pidir,sprintf(...
                    'rfx_-log10pperm_%s.nii',meanres.rows_contrast{con}));
                resvol.data2file(-log10(meanres.pperm(con,:)),pperm);
                pfwe = maxstatpfwe(squeeze(nulldist.r(con,:,:))');
                pfwepath = fullfile(pidir,sprintf(...
                    'rfx_-log10pFWEperm_%s.nii',...
                    meanres.rows_contrast{con}));
                resvol.data2file(-log10(pfwe),pfwepath);
            end
        end
        aap = aas_desc_outputs(aap,'pilab_rsa_r_rfx',rpaths);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
