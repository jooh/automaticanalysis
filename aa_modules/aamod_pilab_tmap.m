% Fit T maps for each chunk in the epi / design instances
% [aap,resp]=aamod_pilab_tmap(aap,task,subj)
function [aap,resp]=aamod_pilab_tmap(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the model / epi instances
        designpath = aas_getfiles_bystream(aap,subj,'pilab_design');
        designvol = load(designpath);
        designvol = designvol.designvol;
        epipath = aas_getfiles_bystream(aap,subj,'pilab_epi');
        epivol = load(epipath);
        epivol = epivol.epivol;
        % find correct labels
        labinds = findStrInArray(designvol.desc.features.unique.labels,...
            aap.tasklist.currenttask.settings.targetname)';
        assert(~isempty(labinds),'found no labels matching %s',...
            aap.tasklist.currenttask.settings.targetname);
        % iterate over chunks
        datcell = cell(epivol.desc.samples.nunique.chunks,1);
        for chind = 1:epivol.desc.samples.nunique.chunks
            c = epivol.desc.samples.unique.chunks(chind);
            datcell{chind} = tmapvol(designvol(...
                designvol.meta.samples.chunks==c,...
                designvol.meta.features.chunks==c),epivol(...
                epivol.meta.samples.chunks==c),labinds);
            % set meta data to avoid problems when concatenating
            datcell{chind}.meta.samples.order(1:datcell{chind}.nsamples,1) = c;
            datcell{chind}.meta.samples.chunks(1:datcell{chind}.nsamples,1) = c;
        end
        % finally, make a big vol
        vol = vertcat(datcell{:});
        % save and describe
        outdir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        outpath = fullfile(outdir,'tvol.mat');
        % very likely too big for older Matlab formats
        save(outpath,'vol','-v7');
        aap=aas_desc_outputs(aap,subj,'pilab_volume',outpath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
