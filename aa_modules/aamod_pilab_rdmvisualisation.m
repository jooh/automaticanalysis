% Visualise a disvol.
% [aap,resp]=aamod_pilab_rdmvisualisation(aap,task,subj)
function [aap,resp]=aamod_pilab_rdmvisualisation(aap,task,subj)

resp='';

switch task
    case 'doit'
        % get mean data RDM
        vpath = aas_getfiles_bystream(aap,subj,'pilab_data_rdms_mean');
        disvol = loadbetter(vpath);
        % get stimuli
        spath = aas_getfiles_bystream(aap,subj,'pilab_stimuli');
        stimuli = loadbetter(spath);

        roinames = disvol.meta.features.names;
        nroi = disvol.nfeatures;

        % prepare output dirs
        pidir = fullfile(aas_getsubjpath(aap,subj),'pilab');
        resdir = fullfile(pidir,'results');
        mkdirifneeded(resdir);
        figdir = fullfile(resdir,'figures');
        mkdirifneeded(figdir);

        ts = aap.tasklist.currenttask.settings;
        % quick check to avoid accidentally computing this on a full
        % searchlight disvol
        assert(nroi <= ts.maxn,'number of ROIs exceed maxn (%d)',ts.maxn);

        % make plots for each ROI
        for roi = 1:nroi
            rankrdm = asrdmmat(tiedrank(disvol.data(:,roi)));
            runplots(rankrdm,stimuli,ts,roinames{roi},figdir,...
                'rank-transformed dissimilarity',{'low','high'});
        end

        % analyse second-order dissimilarity RDM across ROIs
        rdmstim = struct('image',roinames,'alpha',[]);
        rdmat_so = squareform(pdist(disvol.data','spearman'));
        % run plots on that
        runplots(rdmat_so,rdmstim,ts,'second-order rdm',figdir,...
            {'dissimilarity','(Spearman rho)'},[]);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
end

function runplots(rankrdm,stimuli,ts,name,figdir,cblabel,ticklabel)

    printname = strrep(name,' ','');

    % save the raw RDM image
    im = intensity2rgb(rankrdm,ts.cmap);
    imwrite(im,fullfile(figdir,sprintf('rdm_raw_%s.png',...
        printname)),'PNG');
    % save the RDM with stimuli
    F = plotrdms(rankrdm,'labels',{stimuli.image},...
        'nrows',ts.nrows,'titles',name,'cmap',ts.cmap);
    printstandard(fullfile(figdir,sprintf('rdm_wlabel_%s',...
        printname)));
    close(F);
    % same but with colorbar
    F = plotrdms(rankrdm,'labels',{stimuli.image},...
        'nrows',ts.nrows,'colorbarargs',{'label',cblabel,'tick',...
        'minmax','ticklabel',ticklabel},'titles',...
        name,'cmap',ts.cmap,'docb',true);
    printstandard(fullfile(figdir,sprintf('rdm_wcb_%s',...
        printname)));
    close(F);
    % show MDS with stimulus labels, make Shepard plot
    [F,Fshep] = plotrdms(rankrdm,'labels',{stimuli.image},...
        'imagealpha',{stimuli.alpha},'nrows',ts.nrows,'domds',true,...
        'mdsshepard',true,'alphacolor',ts.alphacolor,'mdstimsize',...
        ts.mdstimsize,'titles',name,'dordm',false);
    printstandard(fullfile(figdir,sprintf('mds_%s',...
        printname)),F);
    printstandard(fullfile(figdir,sprintf(...
        'diagnostic_mds_shepard_%s',printname)),Fshep);
    close(F);
    close(Fshep);
end
