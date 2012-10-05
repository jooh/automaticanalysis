% Given a cell/char array of files and a prefix, add the prefix to the filename
% of each
% outfiles = addprefix(files,prefix)
function outfiles = addprefix(files,prefix)

[nfiles,rowlen] = size(files);
outfiles = char([nfiles,rowlen+length(prefix)]);
for f = 1:nfiles
    [path, name, ext] = fileparts(files(f,:));
    outfiles(f,:) = fullfile(path,[prefix name ext]);
end
