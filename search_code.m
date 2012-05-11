function filesFound = search_code(snippet, editFile)
if nargin < 2
    editFile = 0;
end

% We need to be inside the toolbox to work on it
cd(fileparts(mfilename('fullpath')))

toolboxPath = pwd;

fldrDir = genpath(toolboxPath);
addpath(fldrDir); % To add the path to this toolbox!

ind = 0;

filesFound = '';
filters = {'*.m' '*.xml'};

% Then recurse inside each directory until you run out of paths
while ~isempty(strtok(fldrDir, ':'))
    % Get each of the directories made by gendir
    [fldrCurr fldrDir] = strtok(fldrDir, ':');
    
    for f = 1:length(filters)
        
        D = dir(fullfile(fldrCurr, filters{f}));
        for d = 1:length(D)
            T = textread(fullfile(fldrCurr, D(d).name), '%s', 'whitespace', '', 'bufsize', 1024^2);
            
            if ~isempty(strfind(T{1}, snippet))
                filesFound = strvcat(filesFound, sprintf('%s\n', fullfile(fldrCurr, D(d).name)));
                if editFile
                    edit(fullfile(fldrCurr, D(d).name))
                end
            end
        end
    end
end
