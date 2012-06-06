% AA module
% [aap,resp]=aamod_tSNR_EPI(aap,task,subj,sess)
% Calculate the temporal SNR of a series of images within a particular ROI
% Also plots the timecourses...
% This is more accurate (but slower!) than the SNR calculation that is done
% during epi conversion, since it removes signal drift!

function [aap,resp]=aamod_tSNR_EPI(aap,task,subj,sess)

resp='';

switch task        
    case 'doit'
        
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
        ROIimg = aas_getfiles_bystream(aap,subj,'rois');
        
        %% Get started with the processing
        % ROIvol{r}based measures
        EPIsnROI= cell(size(ROIimg,1),1);
        EPIsnHist = cell(size(ROIimg,1),1);
        mROI = cell(size(ROIimg,1),1);
        sROI = cell(size(ROIimg,1),1);
        SNmROI= cell(size(ROIimg,1),1);
        ROIname = cell(size(ROIimg,1),1);
        ROIvol = cell(size(ROIimg,1),1);
        
        % Voxel based measures
        V = spm_vol(deblank(EPIimg(1,:))); % A typical volume...
        EPIsignal = zeros(V.dim(1), V.dim(2), V.dim(3));
        EPInoise = zeros(V.dim(1), V.dim(2), V.dim(3));
        
        fprintf('\nWorking on session %s', aap.acq_details.sessions(sess).name)
        
        fprintf('\n\tLoading ROIs')
        ROIlist = {};
        for r = 1:size(ROIimg,1)
            ROIlist = [ROIlist ROIimg(r,:)];
            [~, ROIname{r}] = fileparts(ROIimg(r,:));
            
            % Now load each of the ROIs we wish to examine (usually using the grey matter)
            rV = spm_vol(ROIimg(r,:));
            ROIvol{r} = spm_read_vols(rV);
            ROIvol{r} = round(ROIvol{r});
            % Mean ROI value is a vector depending on number of scans
            mROI{r} = zeros(size(EPIimg,1),1);
            
            if any(size(ROIvol{r})~=size(EPIsignal))
                aas_log(aap, true, ['The dimensions of the EPI data and the ROI do not match\n' ...
                    'This is likely because you are using normalised EPIs, whereas you should use native ones'])
            end
        end
        
        %% If the dataset is too large, we process it by chunks...
        fprintf('\n\tProcessing data (%d scans)', size(EPIimg,1))
        
        taskComplete = 0;
        chunkDim = 1;
        
        while taskComplete == 0
            fprintf('\nTrying with %d chunks', chunkDim)
            
            try
                chunkX = 0;
                chunkY = 0;
                chunkZ = 0;
                for c = 1:chunkDim
                    chunkX = [chunkX floor(V.dim(1) * c / chunkDim)];
                    chunkY = [chunkY floor(V.dim(2) * c / chunkDim)];
                    chunkZ = [chunkZ floor(V.dim(3) * c / chunkDim)];
                end
                
                % Chunking...
                for x = 1:length(chunkX) - 1
                    for y = 1:length(chunkY) - 1
                        for z = 1:length(chunkZ) - 1
                            fprintf('\n\t...chunk %d %d %d', x, y, z)
                            Xind = chunkX(x) + 1 : chunkX(x+1);
                            Yind = chunkY(y) + 1 : chunkY(y+1);
                            Zind = chunkZ(z) + 1 : chunkZ(z+1);
                            
                            EPIdata = zeros(size(EPIimg,1), length(Xind), ...
                                length(Yind), ...
                                length(Zind));
                            
                            % Load each image into 4-D matrix
                            for e = 1:size(EPIimg,1)
                                V = spm_vol(deblank(EPIimg(e,:)));
                                Y = spm_read_vols(V);
                                EPIdata(e,:,:,:) = Y(Xind,Yind,Zind);
                                
                                %% We can do ROIvol{r} processing here...
                                if x == 1 && y == 1 && z == 1
                                    for r = 1:size(ROIimg,1)
                                        tmp = Y(ROIvol{r}>0);
                                        tmp = tmp(tmp>0); % We don't want to include zero values...
                                        mROI{r}(e,1) = mean(tmp(:)); % Mean per time point
                                        sROI{r}(e,1) = std(tmp(:))./sqrt(length(tmp(:))); % Standard error per time point
                                    end
                                end
                            end
                            
                            % Calcultate signal as mean of the data across volumes
                            EPIsignal(Xind,Yind,Zind) = squeeze(mean(EPIdata, 1));
                            % Calculate noise as standard deviation across volumes
                            EPInoise(Xind,Yind,Zind) = squeeze(std(EPIdata, [], 1));
                        end
                    end
                end
                % If we get here, then we completed the task...
                taskComplete = 1;
            catch aa_error
                %disp(tSNR_error)
                
                if chunkDim > 3
                    error('Error is probably not due to MEMORY')
                end
                
                chunkDim = chunkDim + 1;
            end
        end
        
        fprintf('\n\tCalculating & saving the tSNR image')
        % Calculate SNR as ratio of the two...
        EPIsnr = zeros(size(EPIsignal));
        EPIsnr(EPIsignal ~= 0) = EPIsignal(EPIsignal ~= 0) ./ EPInoise(EPIsignal ~= 0);
        EPIsnr(isnan(EPIsnr)|isinf(EPIsnr)) = 0;
        
        % Save the SNR image!
        sV = V;
        sV.fname = fullfile(aas_getsubjpath(aap,subj), ...
            ['tSNR_' aap.acq_details.sessions(sess).name '.nii']);
        spm_write_vol(sV, EPIsnr);
        
        fprintf('\n\tFinalising ROI data')
        
        for r = 1:size(ROIimg,1)
            % Now get the voxels specific to each ROI
            EPIsnROI{r} = EPIsnr(ROIvol{r}>0);
            EPIsnROI{r} = EPIsnROI{r}(EPIsnROI{r}>0); % We don't want to include zero values...
            % Also get a whole ROI signal and noise estimate
            SNmROI{r} = mean(mROI{r}) ./ std(mROI{r});
        end
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        
        %% tSNR results figure!
        fprintf('\nDisplaying the results of the tSNR analysis')
        colorsB = {'r' 'g' 'b' 'c' 'm' 'y' 'w'};
        
        % We need to make a string for eval, that will print the legend...
        legStr = 'legend(';
        for r = 1:size(ROIimg,1)
            legStr = [legStr ...
                'sprintf(''%s; mn=%.2f; SD=%.2f; med=%.0f; ROI=%.2f (%.0fv)'', ' ...
                'ROIname{' num2str(r) '}, '  ...
                'mean(EPIsnROI{' num2str(r) '}), ' ...
                'std(EPIsnROI{' num2str(r) '}), ' ...
                'median(EPIsnROI{' num2str(r) '}), ' ...
                'SNmROI{' num2str(r) '}, ' ...
                'length(EPIsnROI{' num2str(r) '})),'];
        end
        legStr = [legStr(1:end-1) ');'];
        
        try close(2); catch; end
        
        figure(2)
        set(2, 'Position', [0 0 1200 700])
        maxI = 0;
        windI = 0;
        maxV = 0;
        hold on
        
        for r = 1:size(ROIimg,1)
            % What range do the SNR values take?
            maxI = max(max(EPIsnROI{r}), maxI);
            % What window do we wish to present?
            windI = max(median(EPIsnROI{r}) + std(EPIsnROI{r}) * 3, windI);
        end
        vals = 0:windI/100:ceil(maxI);
        for r = 1:size(ROIimg,1)
            % Now make a histogram and "normalise" it
            EPIsnHist{r} = hist(EPIsnROI{r}, vals);
            EPIsnHist{r} = EPIsnHist{r}./sum(EPIsnHist{r});
            % And decide what is the greatest prop value
            maxV = max(max(EPIsnHist{r}), maxV);
            
            % Make bars semi-transparent for overlaying
            B = bar(vals, EPIsnHist{r}, 1, colorsB{r});
            ch = get(B,'child');
            set(ch, 'faceA', 0.3, 'edgeA', 0.2);
        end
        vals = 0:windI/100:ceil(windI);
        % Set the axis to a good value!
        axis([vals(1), vals(end), 0, maxV*1.05])
        xlabel('SNR')
        ylabel('Proportion of voxels')
        set(gca,'XTick', 0:ceil(maxI./50):maxI)
        title(sprintf('\nSNR for session %s, using %.0f scans', ...
            regexprep(aap.acq_details.sessions(sess).name, '[^a-zA-Z0-9]', ''), ...
            size(EPIimg,1)))
        eval(legStr);
        
        set(gcf,'PaperPositionMode','auto')
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_tSNR.jpeg']));
        
        %% Time-course results figure!
        
        fprintf('\nDisplaying the results of the timecourse analysis')
        
        try close(2); catch; end
        figure(2)
        set(2, 'Position', [0 0 1200 700])
        
        legStr = 'legend(';
        
        for r = 1:size(ROIimg,1)
            %subplot(size(ROIimg,1),1,r)
            hold on
            
            % We need to make a string for eval, that will print the legend...
            legStr = [legStr 'sprintf(''%s (%.0fv)'', ' ...
                'ROIname{' num2str(r) '}, ' ...
                'length(EPIsnROI{' num2str(r) '})),'];
                        
            % Plot main results (errorbars displayed differently now...)
            plot(mROI{r}, ['.' colorsB{r}])
        end
        legStr = [legStr(1:end-1) ');'];
        eval(legStr);
        
        for r = 1:size(ROIimg,1)
            plot(mROI{r} - sROI{r}, '--k')
            plot(mROI{r} + sROI{r}, '--k')
        end
        
        xlim([0 size(EPIimg,1)])
        %ylim([mean(mROI{r} - 2*mean(sROI{r})) mean(mROI{r} + 2*mean(sROI{r}))])
        xlabel('Scan')
        ylabel('Mean signal')
        set(gca,'XTick', 0:ceil(size(EPIimg,1)./25):size(EPIimg,1))
        title(sprintf('\nTimecourse for session %s, using %.0f scans', ...
            regexprep(aap.acq_details.sessions(sess).name, '[^a-zA-Z0-9]', ''), ...
            size(EPIimg,1)))
        
        set(gcf,'PaperPositionMode','auto')
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '_timecourse.jpeg']));
        
        %% Diagnostic VIDEO
        if aap.tasklist.currenttask.settings.diagnostic            
            
            aas_image_avi(sV.fname, ...
                ROIlist, ...
                fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' ROIname{r} '.avi']), ...
                2, ... % Axis
                [800 600], ...
                1); % Rotations
            try close(2); catch; end
        end
        
        %% DESCRIBE OUTPUTS
        
        aap=aas_desc_outputs(aap,subj,'tSNR',sV.fname);
end