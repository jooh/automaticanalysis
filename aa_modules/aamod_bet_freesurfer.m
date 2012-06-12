% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural
% [For best functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_bet_freesurfer(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Let us use the native space...
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            fprintf('\tWARNING: Several structurals found, considering: %s\n', Simg)
        end
        
        [pth nme ext]=fileparts(Simg);
        
        % Structural image (after BETting, if we mask...)
        outStruct=fullfile(pth,['bet_' nme ext]);
        
        % Run MRI_Watershed
        fprintf('1st BET pass (recursive) to find optimal centre of gravity and radius\n')
            
        %% @@@@@
        
        
        %% Make the BET BRAIN MASK
        % As made by BET [and slightly different from inskull_mask]
        V = spm_vol(outStruct);
        M = spm_read_vols(V);
        
        % Mask out non-brain
        M = M > 0;
        
        % Then write out actual BET mask
        V.fname = fullfile(pth, ['bet_' nme '_brain_mask' ext]);
        spm_write_vol(V,M);
        % And add the mask to the list...
        outMask = strvcat(V.fname, outMask);
        
        %% MASK the brain
        fprintf('Masking the brain with Brain Mask \n')
        V = spm_vol(Simg);
        Y = spm_read_vols(V);
        % Mask brain
        Y = Y.*M;
        % Write brain
        V.fname = outStruct;
        spm_write_vol(V, Y);
        
        %% DESCRIBE OUTPUTS!
        
        if aap.tasklist.currenttask.settings.maskBrain
            aap=aas_desc_outputs(aap,subj,'structural',outStruct);
        else
           aap=aas_desc_outputs(aap,subj,'structural',Simg); 
        end
        aap=aas_desc_outputs(aap,subj,'BETmask',outMask);
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        
        %% Draw structural image...
        spm_check_registration(Simg)
        
        % Colour the brain extracted bit pink
        spm_orthviews('addcolouredimage',1,outStruct, [0.9 0.4 0.4])
        
        spm_orthviews('reposition', [0 0 0])
        
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        set(gcf,'PaperPositionMode','auto')
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
        
        %% Diagnostic VIDEO of masks
        if aap.tasklist.currenttask.settings.diagnostic
            
            Ydims = {'X', 'Y', 'Z'};
            for d = 1:length(Ydims)
                aas_image_avi( Simg, ...
                    fullfile(pth, ['bet_' nme '_brain_mask' ext]), ...
                    fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
                    d, ... % Axis
                    [800 600], ...
                    2); % Rotations
            end
            try close(2); catch; end
        end
        
        % Clean up...
        if ~aap.tasklist.currenttask.settings.maskBrain
            delete(outStruct);
        end
end
