% AA module
% Use the Anatomical Transformation Toolbox to normalise the structural to
% a template image

function [aap,resp]=aamod_ANTS_struc2epi(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        %% mean EPI image
        mEPIimg = aas_getfiles_bystream(aap,subj,1,'meanepi');
        
        if isempty(mEPIimg)
            aas_log(aap, true, 'Problem finding mean functional image.');
        elseif size(mEPIimg,1) > 1
            aas_log(aap, false, 'Found more than 1 mean functional images, using first.');
        end
        mEPIimg = deblank(mEPIimg(1,:));
        
        %% Get structural
        % [AVG] Modified the way we get the structural, to be more aa4-like
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural,:));
        
        if isempty(Simg)
            aas_log(aap, true, 'Problem finding structural image.');
        elseif size(Simg,1) > 1
            aas_log(aap, false, 'Found more than 1 structural images, using structural %d', ...
                aap.tasklist.currenttask.settings.structural);
        end
        
        % Get structural directory for this subject
        [Spth, Sfn, Sext] = fileparts(Simg);
        [mEPIpth, mEPIfn, mEPIext] = fileparts(mEPIimg);
        % And make EPI mask if we wish
        if aap.tasklist.currenttask.settings.EPImask
            mask_mEPIimg = fullfile(mEPIpth, ['M' mEPIfn mEPIext]);
            
            % Create a binary mask of the mEPIimg...
            fprintf('Making a binary mEPI mask, to guide deformation!\n')
            V = spm_vol(mEPIimg);
            Y = spm_read_vols(V);
            Y = ones(size(Y));
            V.fname = mask_mEPIimg;
            spm_write_vol(V,Y);
        end
        
        %% Reslice EPI to structural space
        % to keep the structural image in its own space, and avoid slicing
        % off non-EPI acquired data)
        
        % Get realignment defaults
        defs = aap.spm.defaults.realign;
        
        % Flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', 1,...     % what images to reslice
            'mean', 0);           % write mean image
        
        % Reslice
        spm_reslice(strvcat(Simg, mEPIimg), resFlags);
        
        mEPIimg = fullfile(mEPIpth, ['r' mEPIfn mEPIext]);
        
        if aap.tasklist.currenttask.settings.EPImask
            spm_reslice(strvcat(Simg, mask_mEPIimg), resFlags);
            
            mask_mEPIimg = fullfile(mEPIpth, ['rM' mEPIfn mEPIext]);
            
            % Ensure the resliced mask is indeed binary...
            V = spm_vol(mask_mEPIimg);
            Y = spm_read_vols(V);
            Y = Y~=0;
            spm_write_vol(V,Y);
        end
        
        %% Use ANTS to normalise them!
        
        % Set the ANTS path
        setenv('ANTSPATH', aap.directory_conventions.ANTSdir)
        ANTSpath = [fullfile(getenv('ANTSPATH'), 'bin', 'ANTS') ' '];
        warpANTSpath = [fullfile(getenv('ANTSPATH'), 'bin', 'WarpImageMultiTransform') ' '];
        
        % What we get out...
        outfiles = '-o ants.nii ';
        
        % Set up maxiterations
        maxiterations = '';
        for m = 1:length(aap.tasklist.currenttask.settings.maxiterations)
            maxiterations = [maxiterations num2str(aap.tasklist.currenttask.settings.maxiterations(m)) 'x'];
        end
        maxiterations = ['-i ' maxiterations(1:end-1) ' '];
        
        % Dimension number (always 3 for structural)
        Ndim = [num2str(3) ' '];
        
        % Regularisation...
        if ~isempty(aap.tasklist.currenttask.settings.regularisation)
            regularisation = aap.tasklist.currenttask.settings.regularisation;
        else
            regularisation = '';
        end
        
        % SyN transformation...
        if ~isempty(aap.tasklist.currenttask.settings.SyN)
            SyN = ['-t SyN[' num2str(aap.tasklist.currenttask.settings.SyN) '] '];
        else
            SyN = '';
        end
        
        metrics = '';
        for m = 1:9
            if isfield(aap.tasklist.currenttask.settings, ['metric' num2str(m)])
                tmpM = aap.tasklist.currenttask.settings.(['metric' num2str(m)]);
                tmpW = num2str(aap.tasklist.currenttask.settings.(['weight' num2str(m)]));
                tmpP = aap.tasklist.currenttask.settings.(['parameters' num2str(m)]);
                if isnumeric(tmpP)
                    tmpP = num2str(tmpP);
                end
                
                metrics = [ metrics ...
                    '-m ' tmpM '[' mEPIimg ',' Simg ',' tmpW ',' tmpP '] '];
            else
                break
            end
        end
        
        % Any extra options?...
        if ~isempty(aap.tasklist.currenttask.settings.SyN)
            extraoptions = aap.tasklist.currenttask.settings.extraoptions;
        else
            extraoptions = '';
        end
        
        if aap.tasklist.currenttask.settings.EPImask
            EPImask = [' -x ' mask_mEPIimg];
        else
            EPImask = '';
        end
        
        ANTS_command = [ ANTSpath Ndim outfiles maxiterations SyN metrics extraoptions EPImask];
        
        cd(Spth)
        
        % Run ANTS
        fprintf('Running ANTS using command:\n')
        fprintf([ANTS_command '\n'])
        [s w] = aas_shell(ANTS_command);
        
        warpANTS_command = [ warpANTSpath Ndim ... % dimension number
            Simg ' ' fullfile(Spth, ['w' Sfn Sext]) ... % moving image & output
            ' -R ' mEPIimg ' '... % reference image
            fullfile(Spth, 'antsWarp.nii')]; % transform
        if exist(fullfile(Spth,'antsAffine.txt'), 'file')
            warpANTS_command = [warpANTS_command ' ' fullfile(Spth,'antsAffine.txt')]; % and affine, if this exists...
        end
        
        [s w] = aas_shell(warpANTS_command);
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,subj,'structural', fullfile(Spth,['w' Sfn Sext]));
        outANTS = strvcat( ...
            fullfile(Spth,'antsWarp.nii'), ...
            fullfile(Spth, 'antsInverseWarp.nii'));
        if exist(fullfile(Spth,'antsAffine.txt'), 'file')
            outANTS = strvcat(outANTS, fullfile(Spth,'antsAffine.txt'));
        end
        aap=aas_desc_outputs(aap,subj,'ANTs', outANTS);
        
        % Diagnostic image?
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        
        %% Draw coregistration of structural to functional & delta structurals
        V = spm_vol(Simg);
        Y = spm_read_vols(V) - spm_read_vols(spm_vol(fullfile(Spth, ['w' Sfn Sext])));
        V.fname = fullfile(Spth, ['d' Sfn Sext]);
        spm_write_vol(V,Y);
        
        spm_check_registration(strvcat( ...
            fullfile(Spth,['w' Sfn Sext]), ...
            fullfile(Spth,['d' Sfn Sext]), ...
            mEPIimg))
        
        spm_orthviews('reposition', [0 0 0])
        
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
        
        %% Diagnostic VIDEO
        if aap.tasklist.currenttask.settings.diagnostic
            
            Ydims = {'X', 'Y', 'Z'};
            for d = 1:length(Ydims)
                aas_image_avi( mEPIimg, ...
                Simg, ...
                fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
                d, ... % Axis
                [800 600], ...
                2); % Rotations
            end
            try close(2); catch; end
            delete(fullfile(mEPIpth, ['r' mEPIfn mEPIext]))
        end
end
