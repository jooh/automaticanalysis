% AA module wrapping...
% WARNING: This version requires you first to run the aamod_mask_fromstruct
% module present in the MVPaa toolbox and the aamod_fsl_BET module
%
%--------------------------------------------------------------------------
% COMP_SIGNAL creates regressors with mean signal intensity values for each
% segmented compartment [WhiteMatter (WM), CerebralSpinalFluid (CSF) and
% Out-of-Brain (OOB)] separately for each image. The resulting three
% regressors are saved in a mat file (structure: [WM CSF OOB]). These
% regressors can be used similarly to head motion regressors. Where the
% later can account for head motion effects, the compartment signal
% regressors can be used to account for global signal noise due to changes
% in the magnetic field over time (due to the movement of a conductive body
% - like an arm or hand - through the magnetic field) or other nuisance
% factors. Compartment signals are preferred over global signals as
% inclusion of the later might induce fake BOLD deactivations and a
% reduction in power. The compartment signals do not contain GreyMatter (so
% also no BOLD response) and therefore do not suffer from these ill
% effects.
%
% When using these regressors, you could cite my HBM abstract:
%   Verhagen L, Grol MJ, Dijkerman HC, Toni I. (2006) Studying visually-
%       guided reach to grasp movements in an MR-environment. Human Brain
%       Mapping.
%
% or cite my research paper which describes the methods superficially:
%   Verhagen L, Dijkerman HC, Grol MJ, Toni I (2008) Perceptuo-motor
%       interactions during prehension movements. J Neurosci 28(18):4726-4735
%
% or wait for the upcoming methods paper (a little patience required):
%   Verhagen L, Grol MJ, Dijkerman HC, Toni I. Studying visually-guided
%       reach to grasp movements in an MR-environment. Manuscript in
%       preparation.
%
% modified from comp_signal by Lennart Verhagen, 2005-01-25
% L.Verhagen@fcdonders.ru.nl
%--------------------------------------------------------------------------

function [aap,resp]=aamod_compSignal(aap,task,subj,sess)

resp='';

switch task
    case 'domain'
        resp='session';
        
    case 'description'
        resp='Get signal from the CSF, WM, GM and OOB compartments';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Get signal from the CSF, WM, GM and OOB compartments %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        % Let us use the native space...
        SMimg = aas_getfiles_bystream(aap,subj,'segmasksStrict');
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
        BETimg = aas_getfiles_bystream(aap,subj,'epiBETmask');
        
        % Sanity checks
        [junk, fn] = fileparts(EPIimg(1,:));
        indx = strfind(fn, aap.directory_conventions.rawdataafterconversionprefix);
        indx = indx(1);
        if strfind(fn(1:indx-1), 'w')
            aas_log(aap, true, ['You should use unnormalised (i.e. native) images to do this analysis.' ...
                '\n\tThis should be run after the aamod_norm_noss and before aamod_norm_write'])
        end
        
        % Now, let's see which order the masks appear in...
        MOstr = aap.tasklist.currenttask.settings.maskOrder;
        MOlist = {};
        while ~isempty(MOstr)
            [tmp, MOstr] = strtok(MOstr,',');
            MOlist = [MOlist tmp];
        end
        
        % Load the segmented masks!
        mGM = []; mWM = []; mCSF = [];
        for m = 1:size(SMimg,1)
            % The ones *not* including string 'rwc' are the native masks
            [junk, fn] = fileparts(SMimg(m,:));
            if isempty(strfind(fn, 'rwc'))
                indx = strfind(fn, 'rc');
                eval(['m' MOlist{str2num(fn(indx + 2))} ' = spm_read_vols(spm_vol(SMimg(m,:)));'])
            end
        end
        
        % Try to load the BET masks
        for m = 1:size(BETimg,1)
            [junk, fn] = fileparts(BETimg(m,:));
            if strfind(fn, 'outskin_mask')
                mOOH = spm_read_vols(spm_vol(BETimg(m,:)));
                mOOH = ~mOOH; % BET MASK IS INCLUSIVE HEAD...
            elseif strfind(fn, 'skull_mask')
                mSkull = spm_read_vols(spm_vol(BETimg(m,:)));
            end
        end
        % If there is no outskin mask, then OOH will be "not brain"
        if ~exist('mOOH', 'var')
            mOOH = mGM + mWM + mCSF;
            mOOH = mOOH > 0;
            % Mask becomes negative...
            mOOH = ~mOOH;
            fprintf('\nRemoving OOH voxels near Cerebro-Spinal Fluid')
            mOOH = rmNearVox(mOOH, mCSF, aap.tasklist.currenttask.settings.C2Odist);
        end
        
        % Record the number of voxels in each compartment
        nG = sum(mGM(:)>0);
        nW = sum(mWM(:)>0);
        nC = sum(mCSF(:)>0);
        nO = sum(mOOH(:)>0);
        
        if ~isempty(aap.tasklist.currenttask.settings.W2Gdist)
            fprintf('\nRemoving White Matter voxels near Gray Matter')
            mWM = rmNearVox(mWM, mGM, aap.tasklist.currenttask.settings.W2Gdist);
        end
        
        if ~isempty(aap.tasklist.currenttask.settings.C2Gdist)
            fprintf('\nRemoving Cerebro-Spinal Fluid voxels near Gray Matter')
            mCSF = rmNearVox(mCSF, mGM, aap.tasklist.currenttask.settings.C2Gdist);
        end
        
        if ~isempty(aap.tasklist.currenttask.settings.C2Sdist) && exist('mSkull', 'var')
            fprintf('\nRemoving Cerebro-Spinal Fluid voxels near Skull')
            mCSF = rmNearVox(mCSF, mSkull, aap.tasklist.currenttask.settings.C2Sdist);
        end
        
        if ~isempty(aap.tasklist.currenttask.settings.C2Odist)
            fprintf('\nRemoving Cerebro-Spinal Fluid voxels near OOH')
            mCSF = rmNearVox(mCSF, mOOH, aap.tasklist.currenttask.settings.C2Odist);
        end
        
        %% Print the number of voxels in each compartment
        fprintf('\nGrey Matter mask comprises %d (%d) voxels', sum(mGM(:)>0), nG)
        fprintf('\nWhite Matter mask comprises %d (%d) voxels', sum(mWM(:)>0), nW)
        fprintf('\nCereberoSpinal Fluid mask comprises %d (%d) voxels', sum(mCSF(:)>0), nC)
        fprintf('\nOut of Head mask comprises %d (%d) voxels', sum(mOOH(:)>0), nO)
        
        compTC = zeros(size(EPIimg,1), 4);
        for e = 1:size(EPIimg,1)
            Y = spm_read_vols(spm_vol(EPIimg(e,:)));
            % Now average the data from each compartment
            compTC(e,1) = mean(Y(mGM>0));
            compTC(e,2) = mean(Y(mWM>0));
            compTC(e,3) = mean(Y(mCSF>0));
            compTC(e,4) = mean(Y(mOOH>0));
        end
        if any(isnan(compTC(:)))
            aas_log(aap,true, 'Compartment signal contains NaNs')
        end
        
        Rnames = {'GM', 'WM', 'CSF', 'OOH'};
        % Show an image of correlated timecourses...
        corrTCs(compTC, Rnames);
        
        %% DIAGNOSTIC IMAGE
        mriname = aas_prepare_diagnostic(aap,subj);
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
        
        %% Diagnostic VIDEO of masks
        if aap.tasklist.currenttask.settings.diagnostic && ...
                sess == aap.acq_details.selected_sessions(end)
            
            % Write the masks...
            V = spm_vol(BETimg(1,:));
            outSeg = '';
            % GM, WM, CSF, OOH, [Skull]
            V.fname = fullfile(fileparts(BETimg(1,:)), 'GM.nii');
            outSeg = strvcat(outSeg, V.fname);
            spm_write_vol(V, mGM);
            V.fname = fullfile(fileparts(BETimg(1,:)), 'WM.nii');
            outSeg = strvcat(outSeg, V.fname);
            spm_write_vol(V, mWM);
            V.fname = fullfile(fileparts(BETimg(1,:)), 'CSF.nii');
            outSeg = strvcat(outSeg, V.fname);
            spm_write_vol(V, mCSF);
            V.fname = fullfile(fileparts(BETimg(1,:)), 'OOH.nii');
            outSeg = strvcat(outSeg, V.fname);
            spm_write_vol(V, mOOH);
            if exist('mSkull', 'var')
                V.fname = fullfile(fileparts(BETimg(1,:)), 'Skull.nii');
                outSeg = strvcat(outSeg, V.fname);
                spm_write_vol(V, mSkull);
            end
            
            Ydims = {'X', 'Y', 'Z'};
            
            for d = 1:length(Ydims)
                aas_image_avi( BETimg(1,:), ...
                    outSeg, ...
                    fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
                    d, ... % Axis
                    [800 600], ...
                    2, ... % Rotations
                    'fill'); % No outline
            end
            try close(2); catch; end
        end
        
        %% DESCRIBE OUTPUTS!
        
        EPIdir = fileparts(EPIimg(1,:));
        save(fullfile(EPIdir, 'compSignal.mat'), 'compTC')
        aap=aas_desc_outputs(aap,subj,sess,'compSignal',fullfile(EPIdir, 'compSignal.mat'));
end
