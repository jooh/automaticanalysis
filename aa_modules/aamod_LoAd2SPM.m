% AA module
% Use the Anatomical Transformation Toolbox to normalise the structural to
% a template image
% It also outputs a precise BET mask, by using all the segmented tissues

function [aap,resp]=aamod_LoAd2SPM(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        %% Get structural & mask
        % [AVG] Modified the way we get the structural, to be more aa4-like
        SEGimg = aas_getfiles_bystream(aap,subj,'segmentation');
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        
        [Spth, Sfn, Sext] = fileparts(Simg);
        
        tissues = {'WM' 'GM' 'CSF' 'DeepGM' 'Ventricles'};
        tissuesV = cell(size(tissues));
        tissuesY = cell(size(tissues));
        
        for s = 1:size(SEGimg,1)
            for t = 1:length(tissues)
                if ~isempty(strfind(deblank(SEGimg(s,:)), tissues{t}))
                    tissuesV{t} =  spm_vol(deblank(SEGimg(s,:)));
                    tissuesY{t} = spm_read_vols(tissuesV{t}); 
                    break
                end
            end
        end
        
        %% Let's do each of the segmentations...
        V = tissuesV{1};
        outSeg = '';
        
        % GM mask
        Y = 0;
        for t = aap.tasklist.currenttask.settings.includeGM
            Y = Y + tissuesY{t};
        end
        V.fname = fullfile(Spth, ['c1' Sfn Sext]);
        spm_write_vol(V,Y);
        outSeg = strvcat(outSeg, V.fname);
        
        % WM mask
        Y = 0;
        for t = aap.tasklist.currenttask.settings.includeWM
            Y = Y + tissuesY{t};
        end
        V.fname = fullfile(Spth, ['c2' Sfn Sext]);
        spm_write_vol(V,Y);
        outSeg = strvcat(outSeg, V.fname);
             
        % CSF mask
        Y = 0;
        for t = aap.tasklist.currenttask.settings.includeCSF
            Y = Y + tissuesY{t};
        end
        V.fname = fullfile(Spth, ['c3' Sfn Sext]);
        spm_write_vol(V,Y);
        outSeg = strvcat(outSeg, V.fname);
        
        aap=aas_desc_outputs(aap,subj,'segmentation',outSeg);
        
end
end