function aap=aas_checkspmrunning(aap)
% Make sure we get the SPM we really want...
addpath(genpath(aap.directory_conventions.SPMdir))
try
    if (isempty(spm_figure('FindWin')))
        spm('fmri');
    end;
catch
    spm('fmri');
end;
end