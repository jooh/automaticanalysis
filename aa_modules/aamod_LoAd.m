% AA module
% Use the Anatomical Transformation Toolbox to normalise the structural to
% a template image

function [aap,resp]=aamod_LoAd(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        % Set LoAd paths...
        pth=getenv('PATH');
        setenv('PATH',[pth ...
            ':' fullfile(aap.directory_conventions.niftysegdir, 'bin') ...
            ':' fullfile(aap.directory_conventions.niftyregdir, 'bin')]);
        
        %% Get structural & mask
        % [AVG] Modified the way we get the structural, to be more aa4-like
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        BETmask=aas_getfiles_bystream(aap,subj,'BETmask');
        
        [Spth, Sfn, Sext] = fileparts(Simg);
        
        for b = 1:size(BETmask,1)
            if ~isempty(strfind(deblank(BETmask(b,:)), 'brain_mask.nii'))
                BETmask = deblank(BETmask(b,:));
                break
            elseif b == size(BETmask,1)
                aas_log(aap, true, 'Cannot find the correct mask')
            end
        end
        
        %% Use LoAd to segment the structural!
        LoAd_command = ['sh LoAd_brainonly.sh ' ... % Run LoAd command
            Simg ' ' ... % structural
            BETmask]; % mask
        
        [s w] = aas_shell(LoAd_command);
        
        %% Use seg_maths to extract the relevant
        outSeg = '';
        
        tissues = {'WM' 'GM' 'CSF' 'DeepGM' 'Ventricles'};
        
        for t = 1:length(tissues)
            % For each tissue...
            Mfn = fullfile(Spth, [tissues{t} Sext]);
            outSeg = strvcat(outSeg, Mfn);
            segmaths_command = ['seg_maths ' ... % Segment the end result...
                fullfile(Spth, [Sfn '_segmentation' Sext]) ' ' ... % Segmentation image
                '-tp ' num2str(t-1) ' ' ...
                Mfn];
            [s w] = aas_shell(segmaths_command);
        end
        
        aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
        
end
end