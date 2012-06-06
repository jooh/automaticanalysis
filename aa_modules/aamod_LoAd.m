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
        
        
        %% BET mask
        mY = 0;
        % There are 5 sensible tissue classes, the rest are not
        for t = 1:5 
            mY = mY + spm_read_vols(spm_vol(deblank(SEGimg(t,:))));
        end
        mY = mY > 0;
        
        %% @@
        % CAN WE MAKE SURE THAT THERE ARE NO NOISY ISOLATED VOXELS IN BRAIN
        % MASK? TAKE LARGEST CLUSTER...
        %% @@
        
        V.fname = fullfile(Spth, [Sfn '_brain_mask' Sext ]);
        spm_write_vol(V,Y);
        
        % Mask current structural with the more accurate mask, to improve
        % future normalisation
        sV = spm_vol(Simg);
        sY = spm_read_vols(sV);
        sY = sY .* mY;
        spm_write_vol(sV,sY);
        
        %% @@
        % WE NEED SOME NICE DIAGNOSTICS, PERHAPS USING OUR aas_image_avi
        %% @@
        
        % Now put our BETmask in the BETmask stream, but without deleting
        % the original BET mask...
        aap=aas_desc_outputs(aap,subj,'BETmask', ...
            strvcat(V.fname, aas_getfiles_bystream(aap,subj,'BETmask')));
        aap=aas_desc_outputs(aap,subj,'structural',Simg);
        aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
end
end