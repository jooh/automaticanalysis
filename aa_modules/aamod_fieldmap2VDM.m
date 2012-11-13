% AA module
% Convert the fieldmap images (2 mag and 1 phase) into a Voxel Displacement
% Map (VDM) using the FieldMap toolbox of SPM

function [aap,resp]=aamod_fieldmap2VDM(aap,task,subj)

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
        
        if isempty(aap.tasklist.currenttask.settings.tert)
           aas_log(aap, true, 'You must specify a Total EPI Readout Time [in ms]') 
        end
        
        if isempty(aap.tasklist.currenttask.settings.kdir)
            aas_log(aap, true, 'You must specify a Blip Direction [1 or -1] (dependent on PE direction)') 
        end
        
        if isempty(aap.tasklist.currenttask.settings.te1) || ...
                isempty(aap.tasklist.currenttask.settings.te2)
            aas_log(aap, false, 'TE not specified, so let us get it from the fieldmap headers');
           
            HDRfn = aas_getfiles_bystream(aap,subj,'fieldmap_dicom_header');
            
            TE = [];
            HDR = load(HDRfn);
            for h = 1:length(HDR.dcmhdr)
                TE = [TE HDR.dcmhdr{h}.EchoTime];
            end
            TE = unique(TE);
            if length(TE) > 2
                aas_log(aap, true, 'Too many TEs!');
            end
            
            aap.tasklist.currenttask.settings.te1 = TE(1);
            aap.tasklist.currenttask.settings.te2 = TE(2);
        end
        
        
        
        % Defaults specified in this path
        % You can set your own settings in your own copy of the XML or recipe!
        pm_defs = ...
            [aap.tasklist.currenttask.settings.te1, ... % te1
            aap.tasklist.currenttask.settings.te2, ... % te2
            aap.tasklist.currenttask.settings.epifm, ... % epifm
            aap.tasklist.currenttask.settings.tert, ... % tert
            aap.tasklist.currenttask.settings.kdir, ... % kdir
            aap.tasklist.currenttask.settings.mask, ... % (mask)
            aap.tasklist.currenttask.settings.match, ... % (match)
            aap.tasklist.currenttask.settings.writeunwarpedEPI]; % (writeunwarpedEPI)
        
        aas_log(aap, false, sprintf(['Parameters used:' ...
            '\nTE1: %0.3f\tTE2: %0.3f\tTotEPIread: %0.3f\tBlipDir: %d'], ...
            pm_defs(1), pm_defs(2), pm_defs(4), pm_defs(5)));
        
        % Fieldmap path
        FMdir = fullfile(aas_getsubjpath(aap, subj), aap.directory_conventions.fieldmapsdirname);
        
        % The folder fieldmaps must exist...
        if ~exist(FMdir, 'dir')
            mkdir(FMdir)
        end
        
        FMfn = aas_getfiles_bystream(aap,subj,'fieldmap');
        
        % This will work on all sessions (even those we have not selected)
        EPIdir = cell(size(aap.acq_details.sessions, 1));
        for s = aap.acq_details.selected_sessions
            % get files from stream
            EPIdir{s} = aas_getsesspath(aap,subj,s);
        end
        
        % If we cannot find any images in the FMdir, move images there...
        FMfns = dir(fullfile(FMdir, '*.nii'));
        if isempty(FMfns)
            for f = 1:size(FMfn,1)
                unix(['mv ' squeeze(FMfn(f,:)) ' ' FMdir])
            end
        end
        
        FieldMap_preprocess(FMdir,EPIdir,...
            pm_defs,...
            'session');
        
        outstream = {};
        
        % Rename VDM files to their correspondent run names
        if length(EPIdir) == 1
            VDMs = dir(fullfile(FMdir, 'vdm*.nii'));
            [junk, fn, ext] = fileparts(VDMs.name);
            unix(['mv ' fullfile(FMdir, VDMs.name) ' ' ...
                fullfile(FMdir, [fn '_session1' ext])]);
        end
        
        VDMs = dir(fullfile(FMdir, '*session*.nii'));
        
        for v = 1:length(VDMs)
            indx = strfind(VDMs(v).name, 'session');
            s = VDMs(v).name(indx+7:end-4); % Get number after 'session'
            s = str2double(s);
            
            % This gets the selected sessions!
            newfn = [VDMs(v).name(1:indx-1), ...
                aap.acq_details.sessions(aap.acq_details.selected_sessions(s)).name,'.nii'];
            unix(['mv ' fullfile(FMdir, VDMs(v).name)...
                ' ' fullfile(FMdir, newfn)]);
            outstream = [outstream fullfile(FMdir, newfn)];
        end
        
        if isempty(outstream)
            aas_log(aap, true, 'Could not find a fieldmap VDM after processing!')
        end
        
        aap=aas_desc_outputs(aap,subj,'fieldmap',outstream);
        
end
