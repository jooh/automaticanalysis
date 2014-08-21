% visualise RSA RFX results.
%
% [aap,resp]=aamod_pilab_rsa_visualisation_rfx(aap,task)
function [aap,resp]=aamod_pilab_rsa_visualisation_rfx(aap,task)

resp='';

switch task
    case 'doit'
        % get the results
        meanres = loadbetter(aas_getfiles_bystream(aap,'pilab_rsa_r_rfx'));
        groupres = loadbetter(aas_getfiles_bystream(aap,...
            'pilab_rsa_r_group'));
        ts = aap.tasklist.currenttask.settings;

        % prepare output
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        figdir = fullfile(pidir,'figures');

        if ~isempty(ts.pluginpath)
            % call on mean
            feval(ts.pluginpath,figdir,meanres,groupres,...
                'mtarget',ts.mtarget,...
                'errtarget',ts.errtarget,'ptarget',ts.ptarget,...
                'mlabel',ts.mlabel,'errlabel',ts.errlabel,...
                'pthresh',ts.pthresh);
        end

        % standard plots
        plot_roidata(figdir,meanres,groupres,'mtarget',ts.mtarget,...
            'errtarget',ts.errtarget,'ptarget',ts.ptarget,'mlabel',...
            ts.mlabel,'errlabel',ts.errlabel,'pthresh',ts.pthresh);

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
