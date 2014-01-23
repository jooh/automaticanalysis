function [aap,resp]=aamod_sessalign(aap,task,subj)

resp='';

switch task
    case 'doit'

        % figure out the session split
        ts = aap.tasklist.currenttask.settings;
        if isempty(ts.split)
            % just one global split
            ts.split = ones(1,nchunk);
        elseif ischar(ts.split)
            ts.split = eval(ts.split);
        end
        usplit = unique(ts.split);
        nsplit = length(usplit);

        % find the volumes for each split
        for s = 1:nsplit
            thissplit = usplit(s);
            splitinds = find(ts.split==thissplit);
            for r = 1:length(splitinds)
                sess = aap.acq_details.selected_sessions(splitinds(r));
                sessfiles{s}{r} = ...
                    aas_getimages_bystream(aap,subj,sess,'epi');
            end
        end

        [outepi,outmean] = sessalign(sessfiles,aas_getsubjpath(aap,subj));
        
        % Describe outputs
        n = 0;
        for sess = 1:length(outepi)
            for ru = 1:length(outepi{sess})
                n = n+1;
                aap = aas_desc_outputs(aap,subj,n,'epi',outepi{sess}{ru});
            end
        end
        aap = aas_desc_outputs(aap,subj,'meanepi',outmean);
        % nb, no realignment_parameter because I suspect this wouldn't be
        % well defined now
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
