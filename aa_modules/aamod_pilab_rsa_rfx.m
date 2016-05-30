% group analysis of rsa results.
%
% [aap,resp]=aamod_pilab_rsa_rfx(aap,task)
function [aap,resp]=aamod_pilab_rsa_rfx(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);
        % save results in main module directory
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        ts = aap.tasklist.currenttask.settings;

        for s = 1:nsub
            subres(s) = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsa_r'));
            nullres{s} = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsa_nulldist'));
            bootres{s} = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_rsa_bootdist'));
        end
        names = {aap.acq_details.subjects.mriname};
        [subres.name] = names{:};

        fprintf('running roidata_rfx with %d subjects \n',nsub);
        args = {'nperm',ts.nperm,...
            'nboot',ts.nboot,'targetfield',ts.targetfield,...
            'transfun',ts.transfun,'contrasts',ts.contrasts,...
            'minn',ts.minn,'customfits',ts.customfits,'subnull',nullres,...
            'subboot',bootres};
        tic;
        [meanres,groupres] = roidata_rfx(subres,args{:});
        fprintf('finished in %s.\n',seconds2str(toc));

        % save and describe
        outpath = fullfile(pidir,'rsa_r_rfx.mat');
        save(outpath,'meanres');
        aap=aas_desc_outputs(aap,'pilab_rsa_r_rfx',outpath);
        outpath = fullfile(pidir,'rsa_r_group.mat');
        save(outpath,'groupres');
        aap=aas_desc_outputs(aap,'pilab_rsa_r_group',outpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
