% AA module - tsdiffana - tool to assess time series variance
% Rhodri Cusack MRC CBU Cambridge Aug 2004
% subj=subject num
% sess=session num
% Original code written by Doris Eckstein and James Rowe
% Improved (and hopefully not broken) by Rhodri Cusack and Karolina
% Moutsopoulou Jun 2008

function [aap,resp]=aamod_listspikes(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        try close(2); catch; end
        figure(2)
        set(2,'Position', [0 0 800 600])
        
        % lists all scans with high deviation from mean, based on timediff.mat file
        % created in tsdiffana.m, and on rp*.txt file, created by
        % spm_realign.
        
        % Load the movement parameters
        rp = aas_movPars(aap,subj, [1 0 0; 0 0 0]);
        nsess = length(aap.acq_details.sessions);
        
        for sess = aap.acq_details.selected_sessions
            % Load up differnces through time as produced by tsdiffana
            tdfn = aas_getimages_bystream(aap,subj,sess,'tsdiffana');
            
            try
                load (tdfn, 'td', 'globals', 'slicediff');
            catch
                aas_log(aap,1,sprintf('%s not found: Please run tsdiffana first',tdfn));
            end
            
            xyzlimit = aap.tasklist.currenttask.settings.xyzlimit;
            rotlimit_radians=aap.tasklist.currenttask.settings.rotlimit_degrees*pi/180;
            
            %% Now find big changes from one image to the next
            %  td= mean (across voxels) of square difference between one volume and the next
            %  globals= mean global value across an image
            
            tm = td/(mean(globals).^2); % RC/KM added .^2 16/6/2008
            if isempty(aap.tasklist.currenttask.settings.tmlimit)
                % Soft limit
                tmlimit = (mean(tm) + 2 * std(tm));
            else
                % Hard limit
                tmlimit = aap.tasklist.currenttask.settings.tmlimit;
            end
            badimages=[false; (tm > tmlimit)];
            spikes=[find(badimages),tm(badimages(2:end)),slicediff(badimages(2:end))];
            
            %% Now find big movements
            moves = [];
            
            % shift to sync with scan number
            rpdiff = [zeros(1,6); diff(rp{sess})];
            absrpdiff = abs(rpdiff);
            
            badmoves=(any(absrpdiff(:,1:3) > xyzlimit,2)) | (any(absrpdiff(:,4:6) > rotlimit_radians,2));
            moves=[find(badmoves),rpdiff(badmoves,:)];
            
            %% DIAGNOSTIC
            
            % Diagnostic
            subplot(nsess,2,sess * 2 - 1)
            hold on
            plot(tm, 'b')
            plot(spikes(:,1), repmat(tmlimit, [1 size(spikes,1)]), 'r*')
            
            subplot(nsess,2,sess * 2)
            hold on
            plot(rpdiff, 'b')
            plot(moves(:,1), repmat(0, [1 size(moves,1)]), 'r*')
            
            %% Save things
            fprintf('Sess %d \t Spikes: %d; Moves: %d\n', sess, size(spikes,1), size(moves,1))
            
            SPfn = fullfile(aas_getsesspath(aap,subj,sess),sprintf('spikesandmoves.mat'));
            save(SPfn, 'spikes', 'moves');
            
            % Save the time differences
            aap = aas_desc_outputs(aap,subj,sess, 'listspikes', SPfn);
        end
        
        %% Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
