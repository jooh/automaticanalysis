% AA module - tsdiffana - tool to assess time series variance
% [aap,resp]=aamod_tsdiffana(aap,task,subj,sess)
% Rhodri Cusack MRC CBU Cambridge Aug 2004

function [aap,resp]=aamod_tsdiffana(aap,task,subj,sess)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session
    case 'whentorun'
        resp='justonce';  % should this be run everytime or justonce?
    case 'description'
        resp='Run tsdiffana';
    case 'summary'
        resp='Check time series variance using tsdiffana\n';
    case 'report'
        aap.report.html=strcat(aap.report.html,'<table><tr><td>');
        aap=aas_report_addimage(aap,fullfile(aas_getsesspath(aap,subj,sess),'diagnostic_aamod_tsdiffana.jpg'));
        aap.report.html=strcat(aap.report.html,'</td></tr></table>');
    case 'doit'
        
        % get the subdirectories in the main directory
        Spth = aas_getsesspath(aap,subj,sess);
        % get files in this directory
        imgs = aas_getimages_bystream(aap,subj,sess,'epi');
        
        tsdiffana(imgs,0);
        outpath = fullfile(Spth,'timediff.mat');

        % Now produce graphical check
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;

        ts = aap.tasklist.currenttask.settings;
        spikes = [];
        spikeout = fullfile(Spth,'spikes.mat');
        spikethreshout = fullfile(Spth,'spikethresh.mat');
        % any to handle string inputs
        if ~any(isinf(ts.spikethreshold))
            % get the axes for the scaled slice-to-slice variance
            subplot(4,1,1);
            if strcmp(ts.spikethreshold,'manual')
                if exist(spikethreshout,'file')
                    spikethresh = load(spikethreshout);
                    spikethresh = spikethresh.spikethresh;
                    fprintf('loaded existing threshold from %s\n',...
                        spikethreshout);
                else
                    fprintf(['Use the mouse to select a scaled '...
                        'variance threshold for spike identification.\n'])
                    [junk,spikethresh] = ginput(1);
                end
                fprintf('spike threshold = %.2f\n',spikethresh);
            end
            % add thresh to plot
            hold on
            plot(xlim,[spikethresh spikethresh],'r:');
            % vector of spikes
            td = load(outpath);
            % scale similar to plot
            td = td.td / mean(td.globals);
            % +1 since this is a difference measure - so if 1:2 is small
            % and 2:3 is big we want to flag 3, not 2.
            spikes = find(td>spikethresh) + 1;
            fprintf('found %d spikes\n',length(spikes));
            save(spikethreshout,'spikethresh');
        end
        print('-djpeg','-r150',fullfile(Spth,'diagnostic_aamod_tsdiffana'));
        save(spikeout,'spikes');
        
        % Save the time differences
        aap = aas_desc_outputs(aap,subj,sess, 'tsdiffana', outpath);
        aap = aas_desc_outputs(aap,subj,sess, 'spikes', spikeout);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
