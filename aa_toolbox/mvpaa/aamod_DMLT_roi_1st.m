% AA module - MVPaa 1st level (ROI based)
%
% Modified for aa4 by Alejandro Vicente-Grabovetsky Feb-2011

function [aap,resp] = aamod_DMLT_roi_1st(aap,task,p)

resp='';

switch task
    case 'doit'
        
        %% PREPARATIONS
        warning off
        
        mriname = strtok(aap.acq_details.subjects(p).mriname, '/');
        fprintf('Working with data from participant %s. \n', mriname)
        
        % Required...
        aap.tasklist.currenttask.settings.mergeSessions = 1;
        
        %% ANALYSIS
        
        % Load the data into a single big structure...
        [aap data] = mvpaa_loadData(aap, p);
        
        % Load the ROIs from which to extract the data
        try
            ROIimg = aas_getfiles_bystream(aap,p,'rois');
        catch
            % Whole brain?
            ROIimg = aas_getfiles_bystream(aap,p,'segmasksStrict');
            ROIimg = ROIimg(1,:); % [DEBUG] @@@@
        end
        
        ROInum = size(ROIimg,1);
                
        ROIname = {};
        
        DMLout = cell(1,ROInum);
        
        % Loop the routine over all ROIs
        for r = 1:ROInum
            % Get the contrasts for this subject...
            DMLT = mvpaa_loadDMLT(aap,p);
            
            [Rpth Rfn Rext] = fileparts(deblank(ROIimg(r,:)));
            ROIname = [ROIname Rfn];
            
            % Extract betas from each: ROI, voxel, condition, subblock, session
            ROI = uint8(spm_read_vols(spm_vol(fullfile(Rpth, [Rfn Rext]))));
            
            % Check that the ROI size is equal to the data size
            if any(size(ROI) ~= size(data{1,1,1}));
                aas_log(aap, true, 'Your ROI size is different from your data size!');
            end
            
            % Trick for non-binary ROIs...
            if length(unique(ROI))>2
                ROI = ROI > 0;
            end
            voxels = sum(ROI(:));
            
            % ROI to linear index...
            ROI = find(ROI);
            
            Betas = mvpaa_extraction(aap, data, ROI);            
            
            fprintf('\t ROI = %s; vox. = %d (%d)\n',Rfn, sum(~isnan(data{1,1,1}(ROI))), voxels)
            
            if isempty(Betas)
                aas_log(aap, false, sprintf('Not enough voxels in ROI, minimum is %i, you have %i', ...
                    aap.tasklist.currenttask.settings.minVoxels, sum(~isnan(data{1,1,1}(ROI)))));
                continue
            end
            
            for c = 1:length(DMLT)
                con = reshape( ...
                    repmat(DMLT(c).vector(:), [1 aap.tasklist.currenttask.settings.blocks]), ...
                    [1, length(DMLT(c).vector(:))*aap.tasklist.currenttask.settings.blocks])';
                
                Betas = reshape(Betas, [size(Betas,1) size(Betas,2)*size(Betas,3)])';
                
                % The crucial line that calls the DMLT object train method
                DMLT(c).object = DMLT(c).object.train(Betas,con);
                
                % Ideally, the DMLT method you use should already include
                % crossvalidation, etc...
            end
            DMLout{r} = DMLT;
        end
        aap.tasklist.currenttask.settings.ROIname = ROIname;
        
        %% DESCRIBE OUTPUTS
        EP = aap.tasklist.currenttask.settings;
        save(fullfile(aas_getsubjpath(aap,p), [mriname '.mat']), ...
            'DMLout', 'EP')
        aap=aas_desc_outputs(aap,p,'MVPaa', fullfile(aas_getsubjpath(aap,p), [mriname '.mat']));
end