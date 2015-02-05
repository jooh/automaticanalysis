% visualise rsa results. 
%
% [aap,resp]=aamod_pilab_rsa_visualisation(aap,task,subj)
function [aap,resp]=aamod_pilab_rsa_visualisation(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the results
        res = loadbetter(aas_getfiles_bystream(aap,subj,...
            'pilab_result'));

        ts = aap.tasklist.currenttask.settings;
        % prepare output
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        figdir = fullfile(pidir,'figures');

        if ~isempty(ts.pluginpath)
            % call on mean
            arglist = {ts.pluginpath,figdir,meanres,[],...
                'mtarget',ts.mtarget,...
                'errtarget',ts.errtarget,'ptarget',ts.ptarget,...
                'mlabel',ts.mlabel,'errlabel',ts.errlabel,...
                'pthresh',ts.pthresh};
            feval(arglist{:});
        end

        % standard plots
        plot_roidata(figdir,meanres,[],'mtarget',ts.mtarget,...
            'errtarget',ts.errtarget,'ptarget',ts.ptarget,'mlabel',...
            ts.mlabel,'errlabel',ts.errlabel,'pthresh',ts.pthresh);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
