% AAS function that defines finds different files (e.g. niftis) at a determined location

function fileList = aas_fileList(location, prefix, postfix)

if nargin < 1
    error('Need a folder to work on')
end
if nargin < 2
    prefix = '';
end
if nargin < 3
    postfix = '.nii';
end

fileList = {};

D = dir(fullfile(location, [prefix '*' postfix]));

for d = 1:length(D)
    fileList = [fileList, ...
        fullfile(location, ...
        D(d).name)];
end
end