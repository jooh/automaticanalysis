% MVPAA Extraction
% Extracts ROI data from aap.tasklist.currenttask.settingsIs

function Betas = mvpaa_extraction(aap, data, indROI)

voxels = sum(~isnan(data{1,1,1}(indROI)));

% Check that it's worth to extract data
if voxels > aap.tasklist.currenttask.settings.minVoxels
    Betas = nan(voxels, ...
        aap.tasklist.currenttask.settings.conditions, ...
        aap.tasklist.currenttask.settings.blocks, ...
        aap.tasklist.currenttask.settings.sessions);
    
    for s=1:aap.tasklist.currenttask.settings.sessions
        for b=1:aap.tasklist.currenttask.settings.blocks
            for c=1:aap.tasklist.currenttask.settings.conditions
                 tmp = data{c,b,s}(indROI);
                 Betas(:,c,b,s) = tmp(~isnan(tmp));
            end
        end
    end
else
    Betas = [];
end