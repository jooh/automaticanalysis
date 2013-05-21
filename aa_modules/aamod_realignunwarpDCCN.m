% AA module - realignment and unwarp
% As done at the DCCN (Donders Centre for Cognitive Neuroscience)
% [aap,resp]=aamod_realignunwarpDCCN(aap,task,subj)
% Realignment using SPM5
% i=subject num
% Based on aamod_realignunwarp by Rhodri Cusack MRC CBU 2004-6
% Alejandro Vicente Grabovetsky Jan-2012

function [aap,resp]=aamod_realignunwarpDCCN(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        %% Set up a jobs file with some advisable defaults for realign/unwarp!
        jobs = {};
        
        % Get the options from the XML!
        jobs{1}.spatial{1}.realignunwarp.eoptions = ...
            aap.tasklist.currenttask.settings.eoptions;
        jobs{1}.spatial{1}.realignunwarp.uweptions = ...
            aap.tasklist.currenttask.settings.uweoptions;
        jobs{1}.spatial{1}.realignunwarp.uwrptions = ...
            aap.tasklist.currenttask.settings.uwroptions;
                
        % Need to place this string inside a cell?
        jobs{1}.spatial{1}.realignunwarp.eoptions.weight = ...
            {jobs{1}.spatial{1}.realignunwarp.eoptions.weight };
        
        %% Get actual data!
        
        for sess = aap.acq_details.selected_sessions
            fprintf('\nGetting EPI images for session %s', aap.acq_details.sessions(sess).name)
            % Get EPIs
            EPIimg = aas_getimages_bystream(aap,subj,sess,'epi');
            if sess == aap.acq_details.selected_sessions(1)
                % Get first image for the diagnostics...
                diagnosticN = EPIimg(1,:);
            end
            
            jobs{1}.spatial{1}.realignunwarp.data(sess).scans = cellstr(EPIimg);
            
            % Try get VDMs
            try
                % first try to find a vdm with the session name in it
                EPIimg   = spm_select('List', ...
                    fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.fieldmapsdirname), ...
                    sprintf('^vdm.*%s.nii$', aap.acq_details.sessions(sess).name));
                
                % if this fails, try to get a vdm with session%d in it
                if isempty(EPIimg)
                    EPIimg   = spm_select('List', ...
                        fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.fieldmapsdirname), ...
                        sprintf('^vdm.*session%d.nii$',sess));
                end
                vdmpath = fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.fieldmapsdirname, EPIimg);
                assert(exist(vdmpath,'file')>0,'no VDM found!');
                jobs{1}.spatial{1}.realignunwarp.data(sess).pmscan = ...
                    cellstr(vdmpath);
                fprintf('\nFound a VDM fieldmap: %s\n',vdmpath);
            catch
                jobs{1}.spatial{1}.realignunwarp.data(sess).pmscan = ...
                    [];
                fprintf('\nWARNING: Failed to find a VDM fieldmap\n')
            end
        end
        
        %% Run the job!
        
        spm_jobman('run',jobs);
        
        %% Describe outputs
        movPars = {};
        for sess = aap.acq_details.selected_sessions
            uwEPIimg = [];
            for k=1:length(jobs{1}.spatial{1}.realignunwarp.data(sess).scans);
                [pth nme ext] = fileparts(jobs{1}.spatial{1}.realignunwarp.data(sess).scans{k});
                uwEPIimg = strvcat(uwEPIimg,fullfile(pth,['u' nme ext]));
            end
            aap = aas_desc_outputs(aap,subj,sess,'epi',uwEPIimg);
            
            % Get the realignment parameters...
            fn=dir(fullfile(pth,'rp_*.txt'));
            outpars = fullfile(pth,fn(1).name); 
            % Add it to the movement pars...
            movPars = [movPars outpars];
            fn=dir(fullfile(pth,'*uw.mat'));
            outpars = strvcat(outpars, fullfile(pth,fn(1).name));
            aap = aas_desc_outputs(aap,subj,sess,'realignment_parameter',outpars);
            
            if sess == aap.acq_details.selected_sessions(1)
                % mean only for first session
                fn=dir(fullfile(pth,'mean*.nii'));
                aap = aas_desc_outputs(aap,subj,'meanepi',fullfile(pth,fn(1).name));
                
                % Get first image for the diagnostics...
                diagnosticU = uwEPIimg(1,:);
            end
        end
        
        %% Save graphical output to common diagnostics directory
        mriname = aas_prepare_diagnostic(aap,subj);

        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
        
        aas_realign_graph(movPars)
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_MP.jpeg']));
        
        % Let's compare the Native and Unwarped EPI images...
        spm_check_registration(char({diagnosticN; ...
            diagnosticU}))
        
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_EPI.jpeg']));
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
