% RFX group analysis of single subject discriminants.
%
% [aap,resp]=aamod_pilab_rdms_rfx(aap,task)
function [aap,resp]=aamod_pilab_rdms_rfx(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);
        % save results in main module directory
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        ts = aap.tasklist.currenttask.settings;

        for s = 1:nsub
            disvol = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_data_rdms_mean'));
            subres(s) = struct('t',disvol.data,...
                'rows_contrast',{vec2str(1:disvol.nsamples,...
                'dissimilarity_%04d')},'cols_roi',...
                {disvol.meta.features.names},'nfeatures',...
                disvol.meta.features.nfeatures);
        end
        names = {aap.acq_details.subjects.mriname};
        [subres.name] = names{:};

        fprintf('running roidata_rfx with %d subjects \n',nsub);
        tic;
        meanres = roidata_rfx(subres,'nperm',ts.nperm,'nboot',ts.nboot,...
            'targetfield','t');
        fprintf('finished in %s.\n',seconds2str(toc));

        % save as roidata result with p values
        outpath = fullfile(pidir,'rdms_rfx.mat');
        save(outpath,'meanres');
        aap=aas_desc_outputs(aap,'pilab_rdms_rfx',outpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end