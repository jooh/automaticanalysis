% Prepare the diagnostic folder/name, etc.

function mriname = aas_prepare_diagnostic(aap,subj)

% Save graphical output to common diagnostics directory
if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
    mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
end
mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
set(gcf,'PaperPositionMode','auto')