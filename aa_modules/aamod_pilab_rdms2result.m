% Convert a set of disvols to the subres format required for RFX.
%
% [aap,resp]=aamod_pilab_rdms2result(aap,task,subj)
function [aap,resp]=aamod_pilab_rdms2result(aap,task,subj)

resp='';

switch task
    case 'doit'

        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        ts = aap.tasklist.currenttask.settings;
        disvol = loadbetter(aas_getfiles_bystream(aap,subj,...
            'pilab_data_rdms_mean'));
        labels = mat2strcell([1:disvol.nsamples]','dissimilarity_%04d');
        if isfield(disvol.meta.samples,'labels') && ~isempty(disvol.meta.samples.labels)
            labels = disvol.meta.samples.labels;
        end
        subres = struct('d',disvol.data,...
            'rows_contrast',{labels},'cols_roi',...
            {disvol.meta.features.names},'nfeatures',...
            disvol.meta.features.nfeatures,'tail',...
            {repmat({'right'},[disvol.nsamples,1])},...
            'name',aap.acq_details.subjects(subj).mriname);

        % save as roidata result with p values
        outpath = fullfile(pidir,'pilab_result_rdms.mat');
        save(outpath,'subres');
        aap = aas_desc_outputs(aap,subj,'pilab_result',outpath);
        % also need to save dummy null and bootres
        nullres = struct;
        outpath = fullfile(pidir,'pilab_dummy_nulldist.mat');
        save(outpath,'nullres');
        aap = aas_desc_outputs(aap,subj,'pilab_result_nulldist',outpath);
        bootres = struct;
        outpath = fullfile(pidir,'pilab_dummy_bootdist.mat');
        save(outpath,'bootres');
        aap = aas_desc_outputs(aap,subj,'pilab_result_bootdist',outpath);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
