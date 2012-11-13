% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural (usually)
% [For best functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_fsl_FAST(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        firstROIs =  {...
        10, 'L_Thal'; 49, 'R_Thal'; ...
        11, 'L_Caud';  50, 'R_Caud'; ...
        12, 'L_Puta'; 51, 'R_Puta'; ...
        13, 'L_Pall'; 52, 'R_Pall'; ...
        16, 'BrStem'; ...
        17, 'L_Hipp'; 53, 'R_Hipp'; ...
        18, 'L_Amyg'; 54, 'R_Amyg'; ...
        26, 'L_Accu'; 58, 'R_Accu'};
        
        % Find out what stream we should BET
        inputstream = aap.tasklist.currenttask.inputstreams.stream;
        % And the names of the output streams
        outputstream = aap.tasklist.currenttask.outputstreams.stream;
        % And which are the streams which we output...
        outputstream = outputstream(~[strcmp(inputstream,outputstream)]);
        
        % Let us use the native space...
        Simg = aas_getfiles_bystream(aap,subj,inputstream{:});
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            fprintf('WARNING: Several %s found, considering: \n', inputstream{:})
            for t = 1:length(aap.tasklist.currenttask.settings.structural)
                fprintf('\t%s\n', Simg(t,:))
            end
        end
        
        % Image that we will be using for BET...
        [Spth Sfn Sext] = fileparts(Simg);
        
        % Run BET [-R Using robust setting to improve performance!]
        fprintf('Running FSL FIRST\n')
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('run_first_all %s -i %s -o %s', ...
            aap.tasklist.currenttask.settings.options, ...
            Simg, Simg));
        disp(w);
        
        V = spm_vol(fullfile(Spth, [Sfn '_all_fast_firstseg.nii']));
        Y = spm_read_vols(V);
        
        outSeg = {};
        
        % First, we need to create masks for each structure
        for r = 1:length(firstROIs)
            M = Y==firstROIs{r,1};
            if sum(M(:)) > 0
                V.fname = fullfile(Spth, [Sfn '_' firstROIs{r,2} Sext]);
                outSeg = [outSeg, V.fname];
                spm_write_vol(V,M);
            end
        end
        
        %% DESCRIBE OUTPUTS!
        aap=aas_desc_outputs(aap,subj,'rois',outSeg);
        
        %% Save graphical output to common diagnostics directory
        mriname = aas_prepare_diagnostic(aap,subj);
        
        
        
        % This will only work for 1-7 segmentations
        OVERcolours = aas_colours;
        
        %% Draw native template
        spm_check_registration(Simg)
        % Add segmentations...
        for t = 1:length(outSeg)
            spm_orthviews('addcolouredimage',1,outSeg{t}, OVERcolours{t})
        end
        
        spm_orthviews('reposition', [0 0 0])
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
        
        % Another diagnostic image, looking at how well the segmentation worked...
        Pthresh = 0.95;
        
        ROIdata = roi2hist(Simg, ...
            outSeg, Pthresh);
        
        [h, pv, ci, stats] = ttest2(ROIdata{2}, ROIdata{1});
        
        title(sprintf('GM vs WM... T-val: %0.2f (df = %d)', stats.tstat, stats.df))
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_Hist.jpeg']));
        
        %% Diagnostic VIDEO
        if aap.tasklist.currenttask.settings.diagnostic
            Ydims = {'X', 'Y', 'Z'};
            
            for d = 1:length(Ydims)
                if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
                    aas_image_avi(Simg, ...
                        outSeg, ...
                        fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
                        d, ... % Axis
                        [800 600], ...
                        2, ... % Rotations
                        'none'); % No outline...
                    try close(2); catch; end
                end
            end
        end
end
