% AA module
% First-level model Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack Mar 2006-2011
% Additions by Rik Henson Mar 2011

function [aap,resp]=aamod_firstlevel_contrasts(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 contrasts';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Contrasts %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        cwd=pwd;
        % get the subdirectories in the main directory
        subj_dir = aas_getsubjpath(aap,subj);
        
        % Maintained for backwards compatibility- better now now put
        % module-specific value in
        % aap.directory_conventions.stats_singlesubj
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end;
        anadir = fullfile(subj_dir,[aap.directory_conventions.stats_singlesubj stats_suffix]);
        
        % Now set up contrasts...
        SPM=load(aas_getfiles_bystream(aap,subj,'firstlevel_spm'));
        SPM=SPM.SPM;
        SPM.swd=anadir;

        ts = aap.tasklist.currenttask.settings;

        if ts.useaasessions
            sessnames = {aap.acq_details.sessions.name};
            selected_sessions = aap.acq_details.selected_sessions;
            nsess = length(selected_sessions);
            nsess_all = length(sessnames);
        else
            % just get all sessions based on SPM file
            sessnames = {};
            nsess = length(SPM.Sess);
            selected_sessions = 1:nsess;
            nsess_all = nsess;
        end
        
        [fle subjname ext]=fileparts(subj_dir);
        % Load up contrasts from task settings
        contrasts_set=find(strcmp({ts.contrasts.subject},subjname));
        if (isempty(contrasts_set))
            % Try for wildcard
            contrasts_set=find(strcmp({ts.contrasts.subject},'*'));
            if (isempty(contrasts_set))
                aas_log(aap,true,'Can''t find declaration of what contrasts to use - insert this in a local copy of aamod_firstlevel_contrasts.xml or put into user script');
            end;
        end
        contrasts=ts.contrasts(contrasts_set);

        if ~isempty(ts.contrastplugin)
            % use a function plugin to overwrite current cons
            contrasts.con = feval(ts.contrastplugin,subjname);
        end

        % add contrasts for each task regressor v baseline?
        if contrasts.eachagainstbaseline
            basev = zeros(1,length(SPM.Sess(1).col));
            for c = 1:length(basev)
                newv = basev;
                newv(c) = 1;
                contrasts.con(end+1)= struct('name',sprintf(...
                    '%s-o-baseline',SPM.xX.name{SPM.xX.iC(c)}),...
                    'format','sameforallsessions',...
                    'vector',newv,...
                    'session',[],...
                    'type','T');
            end
        end
        
        % logical vector for run-specific contrasts
        % First the general case - across all runs
        nregr = length(SPM.xX.name);
        nruns = length(SPM.Sess);
        runI = true(1,nregr);
        noregr = zeros(1,nregr);
        % Do we also want run-specific contrasts?
        if contrasts.oneconperrun && (nruns > 1)
            for r = 1:nruns
                % All zeros
                runI(r+1,:) = noregr;
                % Then fill in run regs
                runI(r+1,SPM.Sess(r).col)=1;
            end
        end
        
        ccount = 0;
        % Separately for each run config (only one if ~oneconperrun)
        for r = 1:size(runI,1)
            % For each unique contrast
            for conind=1:length(contrasts.con)
                % For each individual contrast (different if oneconperrun)
                ccount = ccount + 1;
                % Get or make automatic contrast name
                if (isfield(contrasts.con(conind),'name') && ~isempty(contrasts.con(conind).name))
                    finalname = contrasts.con(conind).name;
                else
                    finalname = sprintf('Con%d',conind);
                end
                % may have to change name to reflect run
                if r == 1
                    connames{ccount} = finalname;
                else
                    connames{ccount} = sprintf('%s-run%02d',finalname,r-1);
                end
                % support eval'ed strings to define contrasts (e.g. ones, eye)
                if ischar(contrasts.con(conind).vector)
                    contrasts.con(conind).vector = eval(...
                        contrasts.con(conind).vector);
                end
                % Make contract vector
                switch(contrasts.con(conind).format)
                    case {'singlesession','sameforallsessions'}
                        if (strcmp(contrasts.con(conind).format,'singlesession'))
                            sessforcon=[strcmp(sessnames,contrasts.con(conind).session)];
                        else
                            % [AVG] To make the selected sessions work...
                            sessforcon=zeros(1,nsess_all);
                            for sess=selected_sessions
                                sessforcon(sess) = 1;
                            end
                            %sessforcon=ones(1,length(SPM.Sess));
                        end;
                        convec=[];
                        sessnuminspm=1;
                        for sess=selected_sessions
                            numcolsinthissess=length(SPM.Sess(sessnuminspm).col);
                            if (sessforcon(sess))
                                if (size(contrasts.con(conind).vector,2) > numcolsinthissess)
                                    aas_log(aap,true,sprintf('Number of columns in contrast matrix for session %d is more than number of columns in model for this session - wanted %d columns, got ',sess,numcolsinthissess)); disp(contrasts.con(conind).vector);
                                elseif (size(contrasts.con(conind).vector,2) < numcolsinthissess)
                                    convec = [convec contrasts.con(conind).vector zeros(size(contrasts.con(conind).vector,1),numcolsinthissess-size(contrasts.con(conind).vector,2))];
                                    %aas_log(aap,false,sprintf('Warning: Number of columns in contrast matrix for session %d is less than number of columns in model for this session - wanted %d columns, so padding to ',sess,numcolsinthissess)); disp(convec);
                                else
                                    convec=[convec contrasts.con(conind).vector];
                                end
                            else
                                convec=[convec zeros(size(contrasts.con(conind).vector,1),numcolsinthissess)];
                            end;
                            sessnuminspm=sessnuminspm+1;
                        end;
                    case 'uniquebysession'
                        totnumcolsbarconstants = size(SPM.xX.X,2) - nsess;
                        if (size(contrasts.con(conind).vector,2) > totnumcolsbarconstants)
                            aas_log(aap,true,sprintf('Number of columns in contrast matrix for session %d is more than number of columns in model (bar constants) - wanted %d columns, got ',totnumcolsbarconstants)); disp(contrasts.con(conind).vector);
                        elseif (size(contrasts.con(conind).vector,2) < totnumcolsbarconstants)
                            convec = [contrasts.con(conind).vector zeros(size(contrasts.con(conind).vector,1),totnumcolsbarconstants-size(contrasts.con(conind).vector,2))];
                            if (contrasts.automatic_movesandmeans)
                                convec_out=[];
                                ind=1;
                                sessnuminspm=1;
                                for sess=selected_sessions
                                    numcolsinthissess_withoutmoves=length(SPM.Sess(sessnuminspm).col)-6;
                                    newind=ind+numcolsinthissess_withoutmoves;
                                    convec_out=[convec_out convec(:,ind:(newind-1)) zeros(size(convec,1),6)];
                                    ind=newind;
                                    sessnuminspm=sessnuminspm+1;
                                end;
                                convec=convec_out;
                            end;
                            if (size(convec,2) < totnumcolsbarconstants)
                                %aas_log(aap,false,sprintf('Warning: Number of columns in contrast matrix for ''uniquebysession'' option is less than number columns in model (bar constants) = %d, so padding to ',totnumcolsbarconstants)); disp(convec);
                            end;
                        else
                            convec=contrasts.con(conind).vector;
                        end
                    otherwise
                        aas_log(aap,true,sprintf('Unknown format %s specified for contrast %d',contrasts.con(conind).format,ccount));
                end;
                cons{ccount} = [convec zeros(size(convec,1),nsess)];  % Add final constant terms
                
                % Check not empty
                if (~any(cons{ccount}(:)))
                    aas_log(aap,true,sprintf('Contrast %d has no non-zero values, not permitted.',contrasts_set(ccount)));
                end;
                
                % Allow F tests
                if (isfield(contrasts.con(conind),'type') && isempty(contrasts.con(conind).type))
                    contype{ccount}='T';
                else
                    contype{ccount}=contrasts.con(conind).type;
                end;

                % Zero out run-irrelevant entries
                % support for multi-row F contrasts
                nrows = size(cons{ccount},1);
                inds = repmat(runI(r,:)~=1,[nrows 1]);
                cons{ccount}(inds) = 0;
            end;
        end
        
        % Make the con images
        SPM.xCon =[];
        for cc = 1:size(cons,2)
            % skip empty regressors 
            if all(cons{cc}(:) == 0)
                continue
            end
            if length(SPM.xCon)==0
                SPM.xCon = spm_FcUtil('Set',connames{cc},contype{cc},'c',cons{cc}',SPM.xX.xKXs);
            else
                SPM.xCon(end+1) = spm_FcUtil('Set',connames{cc},contype{cc},'c',cons{cc}',SPM.xX.xKXs);
            end
        end
        spm_contrasts(SPM);
        
        % Describe outputs
        %  updated spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));
        
        %  firstlevel_betas (includes related statistical files)
        filters={'con','spmT','spmF'};
        
        for filterind=1:length(filters)
            allbetas=dir(fullfile(anadir,[filters{filterind} '_*']));
            betafns=[];
            for betaind=1:length(allbetas);
                betafns=strvcat(betafns,fullfile(anadir,allbetas(betaind).name));
            end;
            aap=aas_desc_outputs(aap,subj,['firstlevel_' lower(filters{filterind}) 's'],betafns);
        end;
        cd (cwd);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



