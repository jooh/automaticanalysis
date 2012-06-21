% AA module - extended coregistration of EPI to structural
% Coregistration of structural to mean EPI output by realignment in 3 steps
% 1) Coregister Structural to T1 template
% 2) Coregister mean EPI to EPI template
% 3) Coregister mean EPI to Structural
% 4) Apply transformation matrix of mean EPI to all EPIs

function [aap,resp]=aamod_ANTS_build_template(aap,task)

resp='';

switch task
    case 'doit'
        
        % Make template path
        Tpth = fullfile(aas_getstudypath(aap), 'ANTStemplate');
        if ~exist(Tpth, 'dir')
            mkdir(Tpth)
        else
            unix(['rm -rf ' Tpth])
            mkdir(Tpth)
        end
        
        for subj = 1:length(aap.acq_details.subjects)
            Simg = aas_getfiles_bystream(aap,subj,'structural');
            
            % Fileparts to get extension of file...
            [~, ~, Sext] = fileparts(Simg);
            
            if size(Simg,1) > 1
                aas_log(aap, false, 'Found more than 1 structural images, using structural %d', ...
                    aap.tasklist.currenttask.settings.structural);
            end
            % Copy images to right location
            copyfile(Simg, fullfile(Tpth, sprintf('subj%04d%s', subj, Sext)));
        end
        
        %% Use ANTS to make a template!
       
        % Set the ANTS path
        setenv('ANTSPATH', fullfile(aap.directory_conventions.ANTSdir, 'bin/'))
        ANTSpath = [' sh ' fullfile(getenv('ANTSPATH'), 'buildtemplateparallel.sh') ' '];
        % Add the path with functions to interact with torque (qsub) <-- changed
        %setenv('PATH', [getenv('PATH') ':' fullfile(aap.directory_conventions.fieldtripdir, 'qsub')])
        
        % What we get out...
        outfiles = '-o ANTS ';
        
        % Dimension number (always 3 for structural)
        Ndim = ['-d ' num2str(3) ' '];
        
        options = aap.tasklist.currenttask.settings.extraoptions;
        
        ANTS_command = [ ANTSpath Ndim outfiles options ' *'];
        
        cd(Tpth)
        
        % Run ANTS
        fprintf('Running ANTS using command:\n')
        fprintf([ANTS_command '\n'])
        
        [s w] = aas_shell(ANTS_command);
        disp(w)
        
        %% Describe the outputs
        unix(['gunzip ANTStemplate.nii.gz'])
        aap = aas_desc_outputs(aap,'ANTStemplate', fullfile(Tpth, ['ANTStemplate' Sext]));
        
    case 'checkrequirements'
        aas_log(aap,0,'No need to trim or skull strip structural\n' );
end