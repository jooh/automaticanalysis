% visualise discriminant RFX results.
%
% [aap,resp]=aamod_pilab_decode_lindisc_visualisation_rfx(aap,task)
function [aap,resp]=aamod_pilab_decode_lindisc_visualisation(aap,task)

resp='';

switch task
    case 'doit'
        % get the results
        meanres = loadbetter(aas_getfiles_bystream(aap,...
            'pilab_decoder_t_rfx'));
        ts = aap.tasklist.currenttask.settings;

        % prepare output
        pidir = fullfile(aas_getstudypath(aap),'pilab');
        figdir = fullfile(pidir,'figures');

        plot_lindisc(figdir,meanres,'mtarget',ts.mtarget,'errtarget',...
            ts.errtarget,'ptarget',ts.ptarget,'mlabel',ts.mlabel,...
            'errlabel',ts.errlabel,'pthresh',ts.pthresh);
        if ~isempty(ts.pluginpath)
            % call on mean
            feval(ts.pluginpath,figdir,meanres,'mtarget',ts.mtarget,...
                'errtarget',ts.errtarget,'ptarget',ts.ptarget,...
                'mlabel',ts.mlabel,'errlabel',ts.errlabel,...
                'pthresh',ts.pthresh);
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
