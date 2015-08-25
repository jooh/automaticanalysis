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
        maskpath = aas_getfiles_bystream(aap,'pilab_mask_group');
        V = spm_vol(maskpath);
        mask = spm_read_vols(V) > 0;

        parfor s = 1:nsub
            subres(s) = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsa_r'));
            nullres{s} = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsa_nulldist'));
            bootres{s} = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsa_bootdist'));
        end
        names = {aap.acq_details.subjects.mriname};
        [subres.name] = names{:};

        if ~isempty(ts.varsmoothfwhm)
            varsmoothmask = maskpath;
        else
            varsmoothmask = '';
        end

        fprintf('running roidata_rfx with %d subjects \n',nsub);
        tic;
        % NB, no saving of null / bootdist because we don't have the memory
        meanres = roidata_rfx(subres,'nperm',...
            ts.nperm,'nboot',ts.nboot,...
            'targetfield',ts.targetfield,'transfun',ts.transfun,...
            'assumeregister',true,'varsmoothmask',varsmoothmask,...
            'varsmoothfwhm',ts.varsmoothfwhm,'contrasts',ts.contrasts,...
            'customfits',ts.customfits,'subnull',nullres,'subboot',...
            bootres);
        fprintf('finished in %s.\n',seconds2str(toc));

        % save mats
        save(fullfile(pidir,'meanres.mat'),'meanres');

        % save and describe
        % need to dump out each r and p as a nifti
        % probably easiest to first make an SPMVolume instance and then use
        % its methods
        ncon = length(meanres.rows_contrast);
        rpaths = cell(ncon,1);

        mapstowrite(1).maskfun = @(x)[];
        mapstowrite(1).suffix = '';
        mapstowrite(2).maskfun = @(x)x>0;
        mapstowrite(2).suffix = '_less0';
        mapstowrite(3).maskfun = @(x)x<0;
        mapstowrite(3).suffix = '_gt0';

        datadumper = @(data,outpath)datavec2nifti(data,mask,outpath,V);

        for con = 1:length(meanres.rows_contrast)
            % descriptive: mean/stdev
            rpaths{con} = fullfile(pidir,sprintf('rfx_mean_%s.nii',...
                meanres.rows_contrast{con}));
            datadumper(meanres.mean(con,:),rpaths{con});
            datadumper(meanres.std(con,:),fullfile(pidir,...
                sprintf('rfx_std_%s.nii',meanres.rows_contrast{con})));
            if isfield(meanres,'pseudot')
                pseudotpath = fullfile(pidir,sprintf(...
                    'rfx_pseudot_%s.nii',meanres.rows_contrast{con}));
                datadumper(meanres.pseudot(con,:),pseudotpath);
            end

            if strcmp(meanres.tail{con},'both')
                % need to write out three maps for each result - all p, or
                % each side
                thismap = mapstowrite;
            else
                % just write out once
                thismap = mapstowrite(1);
            end

            % write out sign-specific maps possibly
            for m = 1:length(thismap)
                mapind = thismap(m).maskfun(meanres.mean(con,:));
                % parametric p value (with suffix
                pparaout = fullfile(pidir,sprintf(...
                    'rfx_-log10ppara_%s%s.nii',...
                    meanres.rows_contrast{con},thismap(m).suffix));
                ppara = meanres.ppara(con,:);
                % set values that meet test to NaN
                ppara(mapind) = NaN;
                datadumper(-log10(ppara),pparaout);

                % bonferroni p
                pparafweout = fullfile(pidir,sprintf(...
                    'rfx_-log10pFWEpara_%s%s.nii',...
                    meanres.rows_contrast{con},thismap(m).suffix));
                pfwe = ppara * sum(mask(:));
                % no need to set NaNs since we use ppara again
                datadumper(-log10(pfwe),pparafweout);

                % maybe perm p values too?
                if ts.nperm > 1
                    % uncorrected p
                    ppermpath = fullfile(pidir,sprintf(...
                        'rfx_-log10pperm_%s%s.nii',...
                        meanres.rows_contrast{con},thismap(m).suffix));
                    pperm = meanres.pperm(con,:);
                    pperm(mapind) = NaN;
                    datadumper(-log10(pperm),ppermpath);

                    % max stat p
                    pfwepath = fullfile(pidir,sprintf(...
                        'rfx_-log10pFWEperm_%s%s.nii',...
                        meanres.rows_contrast{con},thismap(m).suffix));
                    pfwe = meanres.pfweperm(con,:);
                    pfwe(mapind) = NaN;
                    datadumper(-log10(pfwe),pfwepath);

                    if ~isempty(ts.varsmoothfwhm)
                        % variance smoothed max stat map map
                        pfwesmpath = fullfile(pidir,sprintf(...
                            'rfx_-log10pFWEpermvarsm_%s%s.nii',...
                            meanres.rows_contrast{con},thismap(m).suffix));
                        pfwesm = meanres.pfwepermvarsm(con,:);
                        pfwesm(mapind) = NaN;
                        datadumper(-log10(pfwesm),pfwesmpath);
                    end
                end
            end
        end
        aap = aas_desc_outputs(aap,'pilab_rsa_r_rfx',rpaths);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
