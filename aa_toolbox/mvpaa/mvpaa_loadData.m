% MVPAA Load Data
% Automatically attempts to load data, based on the model you have...

function [aap data] = mvpaa_loadData(aap, subj)

fprintf('Loading beta images \r')

load(aas_getfiles_bystream(aap,subj,'firstlevel_spm'));

% Factor number in each session (must NOT be different)
factorNum = cell(size(SPM.Sess));
% Nuisance conditions
nuisanceNum = cell(size(SPM.Sess));
% block number
blockNum = cell(size(SPM.Sess));
% Conditions, we don't know how many yet...
conditions = cell(size(SPM.Sess));

for s = 1:length(SPM.Sess)
    %% GET CONDITIONS, ETC
    factorNum{s} = nan(size(SPM.Sess(s).U));
    nuisanceNum{s} = zeros(size(SPM.Sess(s).U));
    blockNum{s} = nan(size(SPM.Sess(s).U));
    conditions{s} = cell(size(SPM.Sess(s).U));
    for c = 1:length(SPM.Sess(s).U)
        % Where is the block string for this condition?
        sub = strfind(SPM.Sess(s).U(c).name{:},'sub');
        
        % If condition does not belong to block, then nuisance
        if isempty(sub)
            nuisanceNum{s}(c) = 1;
        else
            % Set number of factors
            factorNum{s}(c) = length(strfind(SPM.Sess(s).U(c).name{:},'_'));
            blockNum{s}(c) = str2double( SPM.Sess(s).U(c).name{:}(sub+3:end));
        end
        indx = strfind(SPM.Sess(s).U(c).name{:}, '_sub');
        if isempty(indx); indx = length(SPM.Sess(s).U(c).name{:}) + 1; end
        conditions{s}{c} = SPM.Sess(s).U(c).name{:}(1:indx-1);
    end
    %% DEAL WITH block NUM
    % Get rid of nuisance conditions...
    blockNum{s}(isnan(blockNum{s})) = [];
    % Find unique blocks
    blocks = unique(blockNum{s});
    % Check that they are equally represented within the session
    for b = 2:length(blocks)
        % If each block contains a different number of conditions...
        if sum(blockNum{s} == blocks(b)) ~= ...
                sum(blockNum{s} == blocks(b-1))
            error('You have a different number of conditions per block')
        end
    end
    blockNum{s} = length(blocks);
    if s > 1
        if blockNum{s} ~= blockNum{s-1}
            error('Unequal number of blocks across the two sessions')
        end
    end
    
    %% DEAL WITH NUISANCE REGRESSORS
    nuisanceNum{s} = sum(nuisanceNum{s});
    if s > 1
        if nuisanceNum{s} ~= nuisanceNum{s-1}
            error('Unequal number of nuisance regressors across the two sessions')
        end
    end
    
    %% DEAL WITH FACTORNUM
    % Get rid of nuisance conditions...
    factorNum{s}(isnan(factorNum{s})) = [];
    
    if all(factorNum{s} == factorNum{s}(1))
        % If number of conditions is constant...
    else
        error('There is something wrong with your regressors...')
    end
    % Now factorNum becomes the number of real different cells
    % Not anymore the number of factors...
    factorNum{s} = sum(factorNum{s}>0)/blockNum{s};
    if s > 1
        if factorNum{s} ~= factorNum{s-1}
            error('Unequal number of factors across the two sessions')
        end
    end
end

%% SAVE SOME PARAMETERS TO AAP
aap.tasklist.currenttask.settings.conditions = factorNum{1};
aap.tasklist.currenttask.settings.blocks = blockNum{1};
aap.tasklist.currenttask.settings.nuisance = nuisanceNum{1};
aap.tasklist.currenttask.settings.sessions = length(SPM.Sess);
% Get rid of repeated condition names...
tempCond = conditions{1}( ...
    (1:aap.tasklist.currenttask.settings.conditions*aap.tasklist.currenttask.settings.blocks) + ...
    aap.tasklist.currenttask.settings.nuisance);
for t1 = length(tempCond):-1:1
    for t2 = t1-1:-1:1
        if strcmp(tempCond{t1}, tempCond{t2})
            tempCond(t1) = [];
            break
        end
    end
end
aap.tasklist.currenttask.settings.conditionNames = tempCond;

%% PREPARATION BEFORE LOADING DATA

% Define data structure
data = cell(aap.tasklist.currenttask.settings.conditions, ...
    aap.tasklist.currenttask.settings.blocks, ...
    aap.tasklist.currenttask.settings.sessions);

fprintf('\nThis experiment contains \n\t%d conditions\n\t%d blocks\n\t%d sessions', ...
    size(data,1), size(data,2), size(data,3))

% Define multiplier depending on the basis function set used...
if strcmp(aap.tasklist.currenttask.settings.basisF, '_TD')
    mult = 3;
elseif strcmp(aap.tasklist.currenttask.settings.basisF, '_T')
    mult = 2;
else
    mult = 1;
end

% Do we grey/white/CSF matter mask the data?
% Get segmentation masks we wish to use, if any

SEGimg = aas_findstream(aap, 'mask', subj);

if ~isempty(SEGimg)
    if aap.tasklist.currenttask.settings.native
        % Select the Segmentation mask we wish to use...
        for m = 1:size(SEGimg,1)
            if strfind(SEGimg(m,:), sprintf('rc%d', aap.tasklist.currenttask.settings.maskNum))
                Mimg = deblank(SEGimg(m,:));
                break
            end
        end
    else
        % Select the Segmentation mask we wish to use...
        for m = 1:size(SEGimg,1)
            if strfind(SEGimg(m,:), sprintf('rwc%d', aap.tasklist.currenttask.settings.maskNum))
                Mimg = deblank(SEGimg(m,:));
                break
            end
        end
    end
    M = spm_read_vols(spm_vol(Mimg));
    % If mask is exclusive, invert it...
    if aap.tasklist.currenttask.settings.maskInclusive == 0;
        M = ~M;
    end
else
    M = [];
end

%% NOW LET US ACTUALLY LOAD DATA (& MASK it with segmentation mask...)

nuis = 0;

Bimg = aas_findstream(aap,'spmts', subj);
datatype = 'spmts';
if ~isempty(Bimg)
    Bimg = aas_findstream(aap,'betas', subj);
    datatype = 'betas';    
end

Bimg = aas_getfiles_bystream(aap,subj,aap.tasklist.currenttask.inputstreams.stream{Sind});

if strcmp(dataType, 'betas')
    prevSess = 0;
    for s = 1:aap.tasklist.currenttask.settings.sessions
        % Works for movement parameters and spikes!
        if s > 1
            prevSess = prevSess + size(SPM.Sess(s - 1).C.C, 2) + length(SPM.Sess(s - 1).U);
        end
        for c=1:aap.tasklist.currenttask.settings.conditions
            for b=1:aap.tasklist.currenttask.settings.blocks
                
                condStr =  [aap.tasklist.currenttask.settings.conditionNames{c} '_sub' num2str(b)];
                
                imageNum = find(strcmp([SPM.Sess(s).U(:).name], ...
                    condStr));
                imageNum = imageNum + prevSess;
                
                % Sanity check to see if conditions have been correctly labelled, etc.
                if isempty(strfind(SPM.Vbeta(imageNum).descrip, condStr))
                    error(['Something went wrong with the condition labelling' ...
                        '\This is probably not your fault! Contact the developer!'])
                end
                
                % Get either betas or or T values...
                V = spm_vol(deblank(Bimg(imageNum*2,:))); % We want .img, not .hdr
                data{c,b,s}=spm_read_vols(V);
                
                if ~isempty(SEGimg)
                    data{c,b,s}(M==0) = NaN;
                end
            end
        end
    end
elseif strcmp(dataType, 'spmts')
    for s = 1:aap.tasklist.currenttask.settings.sessions
        for c=1:aap.tasklist.currenttask.settings.conditions
            for b=1:aap.tasklist.currenttask.settings.blocks
                
                condStr =  [aap.tasklist.currenttask.settings.conditionNames{c} '_sub' num2str(b)];
                
                imageNum = find(strcmp({SPM.xCon.name}, ...
                    condStr));
                imageNum = imageNum(s);
                
                % Get either betas or or T values...
                V = spm_vol(deblank(Bimg(imageNum*2,:))); % We want .img, not .hdr
                data{c,b,s}=spm_read_vols(V);
                
                % Set anything that is 0 to NaN
                data{c,b,s}(data{c,b,s}==0) = NaN;
                if ~isempty(SEGimg)
                    data{c,b,s}(M==0) = NaN;
                end
            end
        end
    end
end