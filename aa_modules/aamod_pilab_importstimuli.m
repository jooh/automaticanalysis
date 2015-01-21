% Import stimulus images using some function.
% [aap,resp]=aamod_pilab_importstimuli(aap,task,subj)
function [aap,resp]=aamod_pilab_importstimuli(aap,task,subj)

resp = '';

switch task
    case 'report'

    case 'doit'
        % find subject name
        subname = aap.acq_details.subjects(subj).mriname;
        ts = aap.tasklist.currenttask.settings;

        stimuli = feval(ts.stimfun,subname);

        if ~isempty(ts.selectinds)
            if ischar(ts.selectinds)
                ts.selectinds = eval(ts.selectinds);
            end
            stimuli = stimuli(ts.selectinds);
        end

        if ~isempty(ts.resortind)
            if ischar(ts.resortind)
                ts.resortind = feval(ts.resortind,numel(stimuli));
            end
            stimuli = stimuli(ts.resortind);
        end

        % save and describe
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        mkdirifneeded(pidir);
        outpath = fullfile(pidir,'pilab_stimuli.mat');
        save(outpath,'stimuli');
        aap = aas_desc_outputs(aap,subj,'pilab_stimuli',outpath);

        % write out as images
        imdir = fullfile(pidir,'images');
        mkdirifneeded(imdir);
        hasalpha = isfield(stimuli,'alpha');
        for s = 1:numel(stimuli)
            im = stimuli(s).image;
            args = {};
            if hasalpha
                args = {'alpha',stimuli(s).alpha};
            end
            imwrite(stimuli(s).image,...
                fullfile(imdir,sprintf('stimulus_%02d.png',s)),'PNG',...
                args{:});
        end

        % Also make a quick diagnostic figure
        fighand = showimages(stimuli,[],ts.plotdims,[],ts.donumbers);
        figout = fullfile(pidir,...
            'diagnostic_pilab_stimuli.png');
        print(figout,'-dpng','-r900');
        close(fighand);
    case 'checkrequirements'
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
