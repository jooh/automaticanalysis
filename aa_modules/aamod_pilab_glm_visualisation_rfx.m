% visualise GLM RFX results.
%
% [aap,resp]=aamod_pilab_glm_visualisation_rfx(aap,task)
function [aap,resp]=aamod_pilab_glm_visualisation_rfx(aap,task)

resp='';

switch task
    case 'doit'
        % get the results
        meanres = loadbetter(aas_getfiles_bystream(aap,'pilab_result_rfx'));
        groupres = loadbetter(aas_getfiles_bystream(aap,...
            'pilab_result_group'));
        ts = aap.tasklist.currenttask.settings;
        arglist = structfields2varargs(ts.roidataargs);

        % prepare output
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        figdir = fullfile(pidir,'figures');

        if ~isempty(ts.pluginpath)
            feval(ts.pluginpath,figdir,meanres,groupres,arglist{:});
        end

        % standard plots
        if ts.runstandardplots
            handles = roidata2figure(meanres,groupres,arglist{:});
            printbyname([handles.figure],figdir);
            close([handles.figure]);
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
