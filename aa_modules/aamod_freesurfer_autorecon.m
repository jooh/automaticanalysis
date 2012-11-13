function [aap,resp]=aamod_freesurfer_autorecon(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Set subject paths
        subjname = aap.acq_details.subjects(subj).mriname;
        subjpath = aas_getsubjpath(aap,subj);
        
        setenv('SUBJECTS_DIR', fileparts(subjpath))
        setenv('FREESURFER_DIR', aap.directory_conventions.freesurferdir)
        
        %% Try to delete old freesurfer running flags
        if exist(fullfile(subjpath, 'ANAT', subjname, 'scripts', 'IsRunning.lh+rh'), 'file')
            unix(['rm ' fullfile(subjpath, 'ANAT', subjname, 'scripts', 'IsRunning.lh+rh')]);
        end
        
        FScommand = ['recon-all -subjid ' subjname ' ' aap.tasklist.currenttask.settings.extraoptions];
        
        disp(FScommand)
        
        [s w] = aas_runFScommand(aap,FScommand);
        
        if aap.tasklist.currenttask.settings.verbose
            disp(w);
        end
        
        if s==1 %|| ~isempty(strfind(w, 'ERROR'))
            error('Some freesurfer ERROR');
        end
        
        %%  make output stream
        FSAR1Dir = fullfile(aas_getsubjpath(aap, subj)); % freesurfer autorecon1 dir
        outstream = dirrec(FSAR1Dir);
        aap=aas_desc_outputs(aap,subj,'freesurfer',outstream);
end
end