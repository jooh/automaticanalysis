% inputs
% sessions: cell array (1 per session) each containing a cell array of file
% paths for EPIs corresponding to each session
% outdir: char path to where results should go
%
% [sessionsout,meanout] = sessalign(sessions,outdir)
function [sessionsout,meanout] = sessalign(sessions,outdir)
orgdir = pwd;


defs = spm_get_defaults;

nsess = length(sessions);
mkdirifneeded(outdir);
outdir_means = fullfile(outdir,'means');
mkdirifneeded(outdir_means);

% rigid body realign each session separately
for s = 1:nsess
    fprintf('realigning session %d of %d\n',s,nsess);
    outdir_sess = fullfile(outdir,sprintf('session_%02d',s));
    mkdirifneeded(outdir_sess);
    cd(outdir_sess);
    spm_realign(sessions{s});
    try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end
    print('-dpng','-r300',fullfile(outdir_sess,'diagnostic_realign'));
    spm_reslice(sessions{s},struct('mean',true));
    % this will be the name of the mean volume that SPM writes out
    meanfiles_org{s} = addprefix(sessions{s}{1}(1,:),'mean');
    [~,fn,ext] = fileparts(meanfiles_org{s});
    meanfiles_re{s} = fullfile(outdir_means,sprintf('session_%02d_%s%s',...
        s,fn,ext));
    % copy over to mean dir
    [success,msg] = copyfile(meanfiles_org{s},meanfiles_re{s});
    assert(success,'copyfile failed: %s',msg);
end

% realign the means
fprintf('realigning the mean EPIs for each session\n')
cd(outdir_means);
spm_realign(meanfiles_re);
spm_reslice(meanfiles_re,struct('mean',true));
% now we expect the template to be
template = addprefix(meanfiles_re{1},'mean');
V_temp = spm_vol(template);

% problem: the qmean files end up beautifully aligned but the EPIs not so
% much. Why? Kind of looks like the data have been normalised but not
% realigned. So all the transforms are there and look correct but the
% volume is in the exact same place where it started.

% general job parameters
spm('defaults','fmri');
spm_jobman('initcfg');
% estimate normalisation between session mean and template and write this
% transform to all realigned session EPI
for s = 1:nsess
    % compile a cell array of all realigned volumes from all runs belonging
    % to this session
    runcell = {};
    sessionsource = {[meanfiles_org{s} ',1']};
    for r = 1:length(sessions{s})
        thisdata = addprefix(sessions{s}{r},defs.realign.write.prefix);
        % this probably assumes you aren't using 4D nifti
        thisdata = [thisdata repmat(',1',[size(thisdata,1) 1])];
        runcell = [runcell cellstr(thisdata)'];
    end
    % add session mean too
    runcell = [runcell sessionsource];
    % session mean as source
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj(s).source = ...
        sessionsource;
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj(s).wtsrc = '';
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj(s).resample = ...
        runcell;
end
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = ...
    {[template ',1']};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
% so far it looks like you get better results without smoothing
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'none';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
% default is 16. It definitely seems to help if you use more.
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 128;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
% these options try to preserve the header info of the template so we can
% overlay in fslview
% Inf is apparently some weird shorthand for use template.
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = Inf([2 3]);
    %world_bb(V_temp);
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = Inf([1 3]);
    %vox2mm(V_temp);
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'q';
% enable post mortems
save(fullfile(outdir,'matlabbatch_job.mat'),'matlabbatch');
spm_jobman('run',matlabbatch);

% figure out what the output files were
sessionsout = sessions;
for s = 1:nsess
    for r = 1:length(sessions{s})
        sessionsout{s}{r} = addprefix(sessionsout{s}{r},['q' ...
            defs.realign.write.prefix]);
        assert(exist(sessionsout{s}{r}(1,:),'file')~=0,...
            'output not found: %s',sessionsout{s}{r}(1,:));
    end
    % keep track of means so we can make a grand average
    V_mean(s) = spm_vol(addprefix(meanfiles_org{s},'q'));
end

% load them up
xyz_mean = spm_read_vols(V_mean);
% allow nans
xyz = mean(xyz_mean,4);
V_out = V_mean(1);
V_out.fname = fullfile(outdir_means,'meanepi.nii');
spm_write_vol(V_out,xyz);
meanout = V_out.fname;

% back to where we started
cd(orgdir);
