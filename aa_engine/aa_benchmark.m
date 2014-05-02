% Automatic analysis - display benchmark
% After some or all of an analysis is run, this command will dump the time
% taken to perform each of the stages
% Rhodri Cusack MRC CBU Cambridge 2004

function [time, flags] = aa_benchmark(studyroot)

% See if we've actually been given an aap, if so retrieve study root
try
    studyroot=fullfile(studyroot.acq_details.root,studyroot.directory_conventions.analysisid);
catch
end;

if (~exist('studyroot','var'))
    studyroot=pwd;
end;

time = []; flags = {};

% First, load AAP structure
aaploadfn=fullfile(studyroot,'aap_parameters');
if exist(aaploadfn,'file'); % djm: prompt if necessary
    aaploadfn=spm_select(1,'mat','Please select aap_parameters file:',pwd);
end
load(aaploadfn);

[s w]=system('hostname');
aas_log(aap,0,['BENCHMARK TIMINGS IN SECONDS ' datestr(now) ' MACHINE ' deblank(w)]);
aas_log(aap,0,['=============================================================']);

% THE MODULES IN AAP.TASKLIST.STAGES ARE RUN IF A CORRESPONDING DONE_ FLAG
% IS NOT FOUND. ONE IS CREATED AFTER SUCCESSFUL EXECUTION
% Now run stage-by-stage tasks
benchmark=[];
benchmark.total=0;
benchmark.max=0;
benchmark.stagesnotdone=0;
for k=1:length(aap.tasklist.main.module)
    aap=aas_setcurrenttask(aap,k);
    deps=aas_dependencytree_allfromtrunk(aap,aap.tasklist.currenttask.domain);
    for depind=1:length(deps)
        doneflag=aas_doneflag_getpath_bydomain(aap,deps{depind}{1},deps{depind}{2},k);
        benchmark=processflag(benchmark, aap, doneflag);
        time = [time; benchmark];
        flags = strvcat(flags, doneflag);
    end;
end;

fprintf('STAGES NOT COMPLETED: %d  PROCESSING TIME: %f  MAX TIME: %f\n',benchmark.stagesnotdone,benchmark.total,benchmark.max);
return;

function [benchmark]=processflag(benchmark,aap, doneflag)
fid=fopen(doneflag,'r');
if (fid==-1)
    benchmark.stagesnotdone=benchmark.stagesnotdone+1;
    fprintf('Stage %s not completed?\n',doneflag);
else
    ln=fgetl(fid);
    tme=str2num(ln);
    try
        ln=fgetl(fid);
        ln=fgetl(fid);
        tme2=str2double(ln);
        tme2str=sprintf('%f',tme2);
        benchmark.max=max(benchmark.max,tme2);
    catch
        tme2str='-';
    end;
    benchmark.total=benchmark.total+tme;
    fprintf('%f\t%s\tfor stage %s\n',tme,tme2str,doneflag);
end;
return;