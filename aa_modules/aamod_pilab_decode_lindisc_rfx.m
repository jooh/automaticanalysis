% group analysis of discriminant results.
%
% [aap,resp]=aamod_pilab_decode_lindisc_rfx(aap,task)
function [aap,resp]=aamod_pilab_decode_lindisc_rfx(aap,task)

resp='';

switch task
    case 'doit'
        nsub = length(aap.acq_details.subjects);
        % save results in main module directory
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        mkdirifneeded(pidir);
        ts = aap.tasklist.currenttask.settings;

        % TODO - outputmode. For now we assume ROI-based. 
        % just load all the results in one go
        for s = 1:nsub
            subres(s) = loadbetter(aas_getfiles_bystream(aap,s,...
                'pilab_decoder_t_mean'));
        end
        names = {aap.acq_details.subjects.mriname};
        [subres.name] = names{:};

        meanres = roidata_rfx(subres,'nperm',ts.nperm,'nboot',ts.nboot,...
            'targetfield','t');

        % save and describe
        outpath = fullfile(pidir,'decoder_t_rfx.mat');
        save(outpath,'meanres');
        aap=aas_desc_outputs(aap,'pilab_decoder_t_rfx',outpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
