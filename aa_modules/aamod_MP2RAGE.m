% AA module
% [aap,resp]=aamod_MP2RAGE(aap,task,subj)
% Take the Inverse Contrast 2 (IC2) and Flat images of an MP2RAGE and obtain a 
% single, BETted (skull stripped) MP2RAGE structural image...
% The IC2 image is useful for BET
% The Flat image contains the image contrast
% This works similarly to bet_robust +
% It corrects bias fields in the structural images!

function [aap,resp]=aamod_MP2RAGE(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';  % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 align';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        % @@@ FIGURE OUT WAY OF MAKING THIS WORK WITH DEFAULT BET FUNCTIONS?
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        DCMfile = aas_getfiles_bystream(aap,subj,'structural_dicom_header');
        
        % dcmhdr{n}.SeriesDescription
        dcmhdr = [];
        load(DCMfile);
        
        % Find the InvC2 and Flat images...
        for d = 1:length(dcmhdr)
            if strfind(dcmhdr{d}.SeriesDescription, 'invContrast2')
                IC2_img = deblank(Simg(d,:));
            elseif strfind(dcmhdr{d}.SeriesDescription, 'flatImg')
                FI_img = deblank(Simg(d,:));
                %{
                % NOT NEEDED!
            elseif strfind(dcmhdr{d}.SeriesDescription, 'invContrast1')
                IC1_fn = Simg(d,:);
            elseif strfind(dcmhdr{d}.SeriesDescription, 'divImg')
                DI_fn = Simg(d,:);
            elseif strfind(dcmhdr{d}.SeriesDescription, 'T1map')
                T1_fn = Simg(d,:);
                %}
            end
        end
        
        Sdir = fileparts(IC2_img);
        
        %% BIAS-FIELD!
        
        % Defaults for normalisation
        defs =aap.spm.defaults.normalise;
        defs.estimate.weight = '';
        
        % Only write out attenuation corrected image
        writeopts.biascor = 1;
        writeopts.GM  = [0 0 0];
        writeopts.WM  = [0 0 0];
        writeopts.CSF = [0 0 0];
        writeopts.cleanup = 0;
        estopts.regtype='';    % turn off affine:
        
        for n = 1:2
            fprintf('Running pass number %d of bias correction\n', n)
            
            out = spm_preproc(IC2_img,estopts);
            sn   = spm_prep2sn(out);
            spm_preproc_write(sn, writeopts);
            
            [pth, fn, ext] = fileparts(IC2_img);
            IC2_img = fullfile(pth, ['m' fn ext]);
        end
        
        %% Do BET on IC2 image
        [pth nme ext]=fileparts(IC2_img);
        
        outStruct=fullfile(pth,['bet_' nme ext]);
        % Run BET [-R Using robust setting to avoid neck!]
        fprintf('1st BET pass (recursive) to find optimal centre of gravity\n')
        [~, w]=aas_runfslcommand(aap, ...
            sprintf('bet %s %s -f %f -v -R',IC2_img,outStruct, ...
            aap.tasklist.currenttask.settings.bet_f_parameter));
        
        % This outputs last radius from recursive command...
        indxS = strfind(w, 'radius');
        indxS = indxS(end) + 7;
        indxE = strfind(w(indxS:end), ' mm');
        indxE = indxE(1) - 2;
        SRad = w(indxS:indxS+indxE);
        
        % We don't extract the centre of gravity from here, since it needs
        % to be input in voxels... Instead get it from betted image
        Y = spm_read_vols(spm_vol(outStruct));
        Y = Y > 0;
        indY = find(Y);
        [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
        COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
        
        fprintf('\t...calculated c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s\n', ...
            COG(1), COG(2), COG(3), SRad)
        
        fprintf('2nd BET pass extracting brain masks in InvCon2 \n')
        % Run BET [-A Now let's get the brain masks and meshes!!]
        [~, w]=aas_runfslcommand(aap, ...
            sprintf('bet %s %s -f %f -c %0.4f %0.4f %0.4f -r %s -v -A',IC2_img,outStruct, ...
            aap.tasklist.currenttask.settings.bet_f_parameter, COG(1), COG(2), COG(3), SRad)...
            );
        
        % This outputs last radius from recursive command...
        indxS = strfind(w, 'radius');
        indxS = indxS(end) + 7;
        indxE = strfind(w(indxS:end), ' mm');
        indxE = indxE(1) - 2;
        SRad = w(indxS:indxS+indxE);
        
        % We don't extract the centre of gravity from here, since it needs
        % to be input in voxels... Instead get it from betted image
        Y = spm_read_vols(spm_vol(outStruct));
        Y = Y > 0;
        indY = find(Y);
        [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
        COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
        
        fprintf('\t...final c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s\n', ...
            COG(1), COG(2), COG(3), SRad)
        
        %% DIAGNOSTIC IMAGE
        % Get the meshes
        D = dir(fullfile(pth, 'bet*mesh*'));
        outMesh = '';
        for d = 1:length(D)
            outMesh = strvcat(outMesh, fullfile(pth, D(d).name));
        end
        
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        
            %% Draw mean EPI...
            spm_check_registration(FI_img)
            
            % This will only work for 1-7 masks
            OVERcolours = {[1 0 0], [0 1 0], [0 0 1], ...
                [1 1 0], [1 0 1], [0 1 1], [1 1 1]};
            
            indx = 0;
            
            try
                % Add mesh outlines, to see if BET has worked properly!
                for r = 1:size(outMesh,1)
                    if strfind(outMesh(r,:), '.nii')
                        indx = indx + 1;
                        spm_orthviews('addcolouredimage',1,outMesh(r,:), OVERcolours{indx})
                    end
                end
            catch
            end
            
            try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.jpeg']));
        
        
        % Mask out the flat image using the BET masked IC2
        
        V = spm_vol(fullfile(pth,['bet_' nme ext]));
        mY = spm_read_vols(V);
        
        V = spm_vol(FI_img);
        fY = spm_read_vols(V);
        
        % Mask out non-brain in fY and write it again
        mY = mY > 0;
        fY = fY.*mY;
        spm_write_vol(V,fY);
        
        % Then write out actual BET mask
        V.fname = fullfile(pth, ['bet_' nme '_brain_mask' ext]);
        spm_write_vol(V,mY);
                    
        % Get the mask images
        D = dir(fullfile(pth, 'bet*mask*'));
        outMask = '';
        for d = 1:length(D)
            outMask = strvcat(outMask, fullfile(pth, D(d).name));
        end
        
        % Get also the meshes
        D = dir(fullfile(pth, 'bet*mesh*'));
        outMesh = '';
        for d = 1:length(D)
            outMesh = strvcat(outMesh, fullfile(pth, D(d).name));
        end
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        try
            %% Draw mean EPI...
            spm_check_registration(IC2_img)
            
            % This will only work for 1-7 masks
            OVERcolours = {[1 0 0], [0 1 0], [0 0 1], ...
                [1 1 0], [1 0 1], [0 1 1], [1 1 1]};
            
            indx = 0;
            
            % Colour the brain extracted bit pink
            spm_orthviews('addcolouredimage',1,FI_img, [0.9 0.4 0.4])
            % Add mesh outlines, to see if BET has worked properly!
            for r = 1:size(outMesh,1)
                if strfind(outMesh(r,:), '.nii')
                    indx = indx + 1;
                    spm_orthviews('addcolouredimage',1,outMesh(r,:), OVERcolours{indx})
                end
            end
            %% Diagnostic VIDEO of masks
            aas_checkreg_avi(aap, subj, 2)
            
            spm_orthviews('reposition', [0 0 0])
            
            try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.jpeg']));
        catch
        end
        
        %% DESCRIBE OUTPUTS!
        
        % Structural image after BETting
        aap=aas_desc_outputs(aap,subj,'structural', FI_img);
        aap=aas_desc_outputs(aap,subj,'BETmask',outMask);
        aap=aas_desc_outputs(aap,subj,'BETmesh',outMesh);        
end