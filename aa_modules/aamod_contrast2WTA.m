% AA module

function [aap,resp]=aamod_contrast2WTA(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='Detect grid orientation from betas coding the sine and cosine';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Grid orientation %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        % Get the SPM structure
        SPMfn = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        anadir = fileparts(SPMfn);
        SPM=load(SPMfn);
        SPM=SPM.SPM;
        
        % Error degrees of freedom
        df = SPM.xX.erdf;
        
        % Winner Takes All matrix
        WTA = 0;
        % Max T-val
        tvalWTA = 0;
        
        tvalY = cell( size(aap.tasklist.currenttask.settings.contrasts));
        
        % Go through each of the contrasts of interest
        for c = aap.tasklist.currenttask.settings.contrasts
            % And the T values corresponding to the contrast
            tvalV = spm_vol(fullfile(anadir, SPM.xCon(c).Vspm.fname));
            tvalY{c==aap.tasklist.currenttask.settings.contrasts} = spm_read_vols(tvalV);
            
            % Get the maximal t value...
            tvalWTA = max(tvalWTA, squeeze(tvalY{c==aap.tasklist.currenttask.settings.contrasts}));
            
            if strcmp('FDR', aap.tasklist.currenttask.settings.threshType)
                [P_thresh, P_thresh_strict] = FDR(1-tcdf(tvalY{c==aap.tasklist.currenttask.settings.contrasts}(:),df), aap.tasklist.currenttask.settings.threshVal);
                T_thresh = abs(tinv(P_thresh,df));
                if isempty(T_thresh)
                    T_thresh = Inf;
                end
            elseif strcmp('uncorr', aap.tasklist.currenttask.settings.threshType)
                T_thresh = abs(tinv(aap.tasklist.currenttask.settings.threshVal,df));
            end
            
            CtvalY = (tvalY{c==aap.tasklist.currenttask.settings.contrasts} >= T_thresh & tvalY{c==aap.tasklist.currenttask.settings.contrasts} >= tvalWTA);
            fprintf('Contrast %s contains %d active voxels\n', SPM.xCon(c).name, sum(CtvalY(:)))
            
            WTA = WTA + CtvalY .* find(aap.tasklist.currenttask.settings.contrasts==c);
        end
        
        V = tvalV;
        
        V.fname = fullfile(anadir, 'WTA.nii');
        spm_write_vol(V,WTA);
        
        WTAimg = V.fname;
        
        %% Now get the ROI and do the average direction & histogram...
        
        ROIimg = aas_findstream(aap,'rois', subj);
        
        if ~isempty(ROIimg)
            
            % Load the ROIs from which to extract the data
            ROIimg = aas_getfiles_bystream(aap,subj,'rois');
            ROI = cell(size(ROIimg,1));
            Rfn = cell(size(ROIimg,1));
            
            for r = 1:size(ROIimg,1)
                [Rpth Rfn{r} Rext] = fileparts(deblank(ROIimg(r,:)));
                
                % Need to extract betas from normalised brains for each: subject, ROI,
                % voxel, epoch, quadrant, subblock
                ROI{r} = int8(spm_read_vols(spm_vol(fullfile(Rpth, [Rfn{r} Rext]))));
                % Trick for non-binary ROIs...
                if length(unique(ROI{r}))>2
                    ROI{r} = ROI{r} > 0;
                end
                ROI{r} = logical(ROI{r});
                voxels = sum(ROI{r}(:));
                
                fprintf('\t ROI = %s; vox. = %d\n',Rfn{r}, voxels)
                
                rWTA = WTA;
                rWTA(~ROI{r}) = 0;
                
                V.fname = fullfile(anadir, ['WTA_' Rfn{r} '.nii']);
                spm_write_vol(V,rWTA);
                
                WTAimg = strvcat(WTAimg, V.fname);
                
                % And some line plots...
                tvalMat = nan(length(aap.tasklist.currenttask.settings.contrasts), ...
                    voxels);
                
                for c = 1:length(aap.tasklist.currenttask.settings.contrasts)
                    tvalMat(c,:) = tvalY{c}(ROI{r});
                end
            end
        end
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        if ~isempty(ROIimg)
            % Display the pattern of voxel-wise T-values across the contrasts
            figure(2)
            set(2, 'Position', [0 0 1000 800])
            
            for r = 1:size(ROIimg,1)
                h = subplot(size(ROIimg,1), 1, r);
                plot(tvalMat, '-')
                
                set(gca, ...
                    'xtick', 1:length(aap.tasklist.currenttask.settings.contrasts), ...
                    'xticklabel', {SPM.xCon(aap.tasklist.currenttask.settings.contrasts).name})
                
                xlabel('Contrasts')
                ylabel('T-values')
                
                % Display the correlation between the T-values across voxels
                figure(3)
                set(3, 'Position', [0 0 1000, 600])
                subplot(1,3,1)
                tvalSimil = squareform(pdist(tvalMat'));
                imagesc(tvalSimil)
                axis equal off
                colorbar
                title('Voxel dissimilarity across regressors')
                
                subplot(1,3,2)
                indROI = find(ROI{r});
                [xROI yROI zROI] = ind2sub(size(ROI{r}), indROI);
                roiSimil = squareform(pdist([xROI yROI zROI]));
                imagesc(roiSimil);
                axis equal off
                colorbar
                title('Voxel distance within ROI')
                
                
                %% CLUSTERING (move out?!)
                CtvalSimil = tvalSimil;
                uniqueSimil = logical(triu(ones(sum(ROI{r}(:))), 1));
                regs = [];
                distPow = [1:6]; %[(1/3) (1/2) 1 2 3];
                for w = distPow
                    regs = [regs roiSimil(uniqueSimil).^w - mean(roiSimil(uniqueSimil).^w)];
                end
                
                [b, dev, stats] = glmfit(regs, ...
                    CtvalSimil(uniqueSimil));
                
                for w = 1:length(distPow)
                    CtvalSimil = CtvalSimil - ...
                        b(1+w).*(roiSimil.^w - mean(roiSimil(:).^w));
                end
                
                subplot(1,3,3)
                imagesc(CtvalSimil);
                axis equal off
                colorbar
                title('Voxel dissimilarity across regressors (corrected)')
                
                fprintf('%s distance to contrast similarity: %0.3f\n', Rfn{r}, corr(tvalSimil(uniqueSimil), roiSimil(uniqueSimil)))
                fprintf('%s distance to contrast similarity: %0.3f (corrected)\n', Rfn{r}, corr(CtvalSimil(uniqueSimil), roiSimil(uniqueSimil)))
                fprintf('Shared variance between original & corrected: %0.3f\n', corr(CtvalSimil(uniqueSimil), tvalSimil(uniqueSimil)).^2)
                
                set(gcf,'PaperPositionMode','auto')
                print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' mriname '_' Rfn{r} '.jpeg']));
                
                % Cluster analysis
                clusterNum = 100;
                CtvalSimil(logical(eye(size(CtvalSimil)))) = 0;
                CtvalLinkage = linkage(squareform(CtvalSimil), 'single');
                CtvalCluster = cluster(CtvalLinkage, 'cutoff', 2, 'criterion', 'inconsistent');
                
                clusterNum = max(CtvalCluster(:));
                
                clusterROI = zeros(size(ROI{r}));
                clusterROI(~ROI{r}) = nan;
                clusterROI(ROI{r}) = 0;
                for c = 1:clusterNum
                    clusterROI(indROI(CtvalCluster==c)) = c;
                end
                
                figure(4)
                for z = 1:size(clusterROI,3)
                    sliceROI = clusterROI(:,:,z);
                    if nansum(sliceROI(:)) > 0
                        imagescnan(sliceROI)
                        caxis([0 clusterNum])
                        axis equal off
                        colorbar
                        pause(1)
                    end
                end                
            end
            
            figure(2)
            set(gcf,'PaperPositionMode','auto')
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '_timecourse.jpeg']));
        end
        
        %% DESCRIBE OUTPUTS!
        aas_desc_outputs(aap,subj,'WTA', WTAimg);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
