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

        % prepare output
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        figdir = fullfile(pidir,'figures');

        arglist = {figdir,meanres,groupres,[],...
            'mtarget',ts.mtarget,...
            'errtarget',ts.errtarget,'ptarget',ts.ptarget,...
            'mlabel',ts.mlabel,'errlabel',ts.errlabel,...
            'pthresh',ts.pthresh,'extracongroups',ts.extracongroups,...
            'groupmtarget',ts.groupmtarget,'groupptarget',...
            ts.groupptarget};
        if ~isempty(ts.pluginpath)
            feval(ts.pluginpath,arglist{:});
        end

        % standard plots
        if ts.runstandardplots
            plot_roidata(arglist{:});
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
