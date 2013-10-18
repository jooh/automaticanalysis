% visualise discriminant results. 
%
% [aap,resp]=aamod_pilab_decode_lindisc_visualisation(aap,task,subj)
function [aap,resp]=aamod_pilab_decode_lindisc_visualisation(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get the results
        meanpath = aas_getfiles_bystream(aap,subj,'pilab_decoder_t_mean');
        meanres = loadbetter(meanpath);
        sesspath = aas_getfiles_bystream(aap,subj,'pilab_decoder_t_sess');

        ts = aap.tasklist.currenttask.settings;

        % prepare output
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
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

        for s = 1:size(sesspath,1)
            sessres = loadbetter(sesspath(s,:));
            sessfig = fullfile(figdir,sprintf('sess%02d',s));
            plot_lindisc(sessfig,sessres,'mtarget',ts.mtarget,...
                'errtarget',ts.errtarget,'ptarget',ts.ptarget,...
                'mlabel',ts.mlabel,'errlabel',ts.errlabel,...
                'pthresh',ts.pthresh);
            if ~isempty(ts.pluginpath)
                feval(ts.pluginpath,sessfig,sessres,'mtarget',...
                    ts.mtarget,'errtarget',ts.errtarget,'ptarget',...
                    ts.ptarget,'mlabel',ts.mlabel,...
                    'errlabel',ts.errlabel,'pthresh',ts.pthresh);
            end
        end

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
