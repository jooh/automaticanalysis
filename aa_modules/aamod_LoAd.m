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
            ':' fullfile(aap.directory_conventions.niftyregdir, 'bin') ...
            ]);
        
        % Set libraries
        libpth = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH', [libpth ...
            ':' fullfile(aap.directory_conventions.niftysegdir, 'lib') ...
            ':' fullfile(aap.directory_conventions.niftyregdir, 'lib') ...
            ]);
        
        % Set libraries (mac)
        libMpth = getenv('DYLD_LIBRARY_PATH');
        setenv('DYLD_LIBRARY_PATH', [libMpth ...
            ':' fullfile(aap.directory_conventions.niftysegdir, 'lib') ...
            ':' fullfile(aap.directory_conventions.niftyregdir, 'lib') ...
            ]);
        
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
        
        % Let us loop the load algorhythm several times, so as to obtain
        % a single WM and GM mass...
        looseBits = 1;
        
        while looseBits == 1
            
            %% Use LoAd to segment the structural!
            LoAd_command = ['sh LoAd_brainonly.sh ' ... % Run LoAd command
                Simg ' ' ... % structural
                BETmask]; % mask
            
            [s w] = aas_shell(LoAd_command);
            disp(w)
            
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
                disp(w)
            end
            
            %% BET mask
            mY = 0;
            % There are 5 sensible tissue classes, the rest are not
            % Exclude CSF, as this will make the brain mask too large...
            for t = [1 2 4 5]
                mY = mY + spm_read_vols(spm_vol(deblank(SEGimg(t,:))));
            end
            mY = mY > 0;
            
            % Remove loose bits
            CC = bwconncomp(mY);
            if length(CC.PixelIdxList) == 1
                looseBits = 0;
            else
                numPixels = cellfun(@numel,CC.PixelIdxList);
                [biggest,idx] = max(numPixels);
                for b = 1:length(CC.PixelIdxList)
                    if b~=inx
                        mY(CC.PixelIdxList{idx}) = 0;
                    end
                end
            end
            
            % Fill any holes that may be remaining
            try
                mY = imfill(mY,'holes')
            catch aa_error
                aas_log(aap, false, 'Could not fill the holes in brain mask')
            end
            
            BETmask = fullfile(Spth, [Sfn '_brain_mask' Sext ]);
            V.fname = BETmask;
            spm_write_vol(V,Y);
        end
        
        % Mask current structural with the more accurate mask, to improve
        % future normalisation
        sV = spm_vol(Simg);
        sY = spm_read_vols(sV);
        sY = sY .* mY;
        spm_write_vol(sV,sY);
        
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        %% Diagnostic VIDEO
        if aap.tasklist.currenttask.settings.diagnostic
            
            Ydims = {'X', 'Y', 'Z'};
            
            for d = 1:length(Ydims)
                aas_image_avi( Simg, ...
                    outSeg, ...
                    fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
                    d, ... % Axis
                    [800 600], ...
                    2, ... % Rotations
                    'none'); % No outline
            end
            try close(2); catch; end
            delete(fullfile(mEPIpth, ['r' mEPIfn mEPIext]))
        end
        
        % Now put our BETmask in the BETmask stream, but without deleting
        % the original BET mask...
        if aap.tasklist.currenttask.settings.brainMask
            aap=aas_desc_outputs(aap,subj,'BETmask', ...
                strvcat(V.fname, aas_getfiles_bystream(aap,subj,'BETmask')));
        else
            delete(BETmask)
        end
        aap=aas_desc_outputs(aap,subj,'structural',Simg);
        aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
end
end