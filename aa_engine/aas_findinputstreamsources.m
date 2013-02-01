
%% CONNECT DATA PIPELINE BY IDENTIFYING SOURCE OF INPUT FROM OUTPUT STREAMS
%   Works back through dependencies to determine where inputs come from
function [aap]=aas_findinputstreamsources(aap)

% Make empty cell structures for input and output streams
aap.internal.inputstreamsources=cell(length(aap.tasklist.main.module),1);
aap.internal.outputstreamdestinations=cell(length(aap.tasklist.main.module),1);
for k1=1:length(aap.tasklist.main.module)
    aap.internal.inputstreamsources{k1}.stream=[];
    aap.internal.outputstreamdestinations{k1}.stream=[];
end;

% Now go through each module and find its input dependencies
%  then make bi-directional connections in inputstreamsources and
%  outputstreamdestinations
for k1=1:length(aap.tasklist.main.module)
    [stagepath stagename]=fileparts(aap.tasklist.main.module(k1).name);
    index=aap.tasklist.main.module(k1).index;
    
    % Find streams to be loaded remotely
    remotestream=aap.tasklist.main.module(k1).remotestream;
    
    if (isfield(aap.schema.tasksettings.(stagename),'inputstreams'))
        inputstreams=aap.schema.tasksettings.(stagename).inputstreams;
        
        for i=1:length(inputstreams.stream)
            inputstreamname=inputstreams.stream{i};
            ismodified=1;
            if isstruct(inputstreamname)
                if isfield(inputstreamname.ATTRIBUTE,'ismodified')
                    ismodified=inputstreamname.ATTRIBUTE.ismodified;
                end;
                inputstreamname=inputstreamname.CONTENT;
            end;
            findremote=[];
            if ~isempty(remotestream)
                findremote=find(strcmp(inputstreamname,{remotestream.stream}));
            end;
            if ~isempty(findremote)
                stream=[];
                stream.name=inputstreamname;
                stream.sourcenumber=-1;
                stream.sourcestagename=remotestream(findremote).stagetag;
                stream.sourcedomain=remotestream(findremote).sourcedomain;
                stream.depth=[];
                stream.host=remotestream(findremote).host;
                stream.aapfilename=remotestream(findremote).aapfilename;
                stream.ismodified=ismodified;
                if (isempty(aap.internal.inputstreamsources{k1}.stream))
                    aap.internal.inputstreamsources{k1}.stream=stream;
                else
                    aap.internal.inputstreamsources{k1}.stream(end+1)=stream;
                end;
                % [AVG] to make the inputs/outpus more obvious!
                %aas_log(aap,false,sprintf('Stage %s input %s comes from remote host %s stream %s',stagename,stream.name,stream.host,stream.sourcestagename));
                % JC - don't use cprintf directly since this causes crashes
                % if running Matlab from a terminal. aas_log handles this
                % case much better
                aas_log(aap,0,'Stage','text');
                aas_log(aap,0,stagename,[65, 105, 225]/255);
                aas_log(aap,0,'input','text');
                aas_log(aap,0,stream.name,[46, 139, 87]/255);
                aas_log(aap,0,'comes from remote host','text');
                aas_log(aap,0,stream.host,'-text');
                aas_log(aap,0,'stream','text');
                aas_log(aap,0,stream.sourcestagename,'Magenta');
            else
                
                [aap stagethatoutputs mindepth]=aas_searchforoutput(aap,k1,inputstreamname,true,0,inf);
                if isempty(stagethatoutputs)
                    % [AVG] to make the inputs/outpus more obvious!
                    %aas_log(aap,true,sprintf('Stage %s required input %s is not an output of any stage it is dependent on. You might need to add an aas_addinitialstream command or get the stream from a remote source.',stagename,inputstreamname));
                    
                    aas_log(aap,0,'Stage','text');
                    aas_log(aap,0,stagename,[65, 105, 225]/255);
                    aas_log(aap,0,'required input','text');
                    aas_log(aap,0,inputstreamname,[46, 139, 87]/255);
                    aas_log(aap,0,'is not an output of any stage it is dependent on','text')
                    aas_log(aap,true,'You might need to add an aas_addinitialstream command or get the stream from a remote source.');
                else
                    [sourcestagepath sourcestagename]=fileparts(aap.tasklist.main.module(stagethatoutputs).name);
                    sourceindex=aap.tasklist.main.module(stagethatoutputs).index;
                    % [AVG] to make the inputs/outpus more obvious!
                    %aas_log(aap,false,sprintf('Stage %s input %s comes from %s which is %d dependencies prior',stagename,inputstreamname,sourcestagename,mindepth));
                    aas_log(aap,0,'Stage','text');
                    aas_log(aap,0,stagename,[65, 105, 225]/255);
                    aas_log(aap,0,'input','text');
                    aas_log(aap,0,inputstreamname,[46, 139, 87]/255);
                    aas_log(aap,0,'comes from','text');
                    aas_log(aap,0,sourcestagename,'Magenta');
                    aas_log(aap,0,'which is','text');
                    aas_log(aap,0, num2str(mindepth),'-text');
                    aas_log(aap,0,'dependencies prior.\n','text')
                    
                    stream=[];
                    stream.name=inputstreamname;
                    stream.sourcenumber=stagethatoutputs;
                    stream.sourcestagename=sourcestagename;
                    stream.sourcedomain=[];
                    stream.depth=mindepth;
                    stream.host='';
                    stream.aapfilename='';
                    stream.ismodified=ismodified;
                    stream.sourcedomain=aap.schema.tasksettings.(sourcestagename)(sourceindex).ATTRIBUTE.domain;
                    if (isempty(aap.internal.inputstreamsources{k1}.stream))
                        aap.internal.inputstreamsources{k1}.stream=stream;
                    else
                        aap.internal.inputstreamsources{k1}.stream(end+1)=stream;
                    end;
                    
                    stream=[];
                    stream.name=inputstreamname;
                    stream.destnumber=k1;
                    stream.deststagename=stagename;
                    stream.depth=mindepth;
                    stream.destdomain=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.domain;
                    if (isempty(aap.internal.outputstreamdestinations{stagethatoutputs}.stream))
                        aap.internal.outputstreamdestinations{stagethatoutputs}.stream=stream;
                    else
                        aap.internal.outputstreamdestinations{stagethatoutputs}.stream(end+1)=stream;
                    end;
                end;
            end;
        end;
    end;
end;
