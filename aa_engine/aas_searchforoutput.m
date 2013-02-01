% RECURSIVELY SEARCH DEPENDENCIES
%  to see which will have outputted each
%  input required for this stage
%  note that inputs can be affected by dependency map
%  (moved out of aas_findinputstreamsources because subfunctions are evil)
function [aap,stagethatoutputs,mindepth]=aas_searchforoutput(aap,currentstage,outputtype,notthislevelplease,depth,mindepth)

% is this branch ever going to do better than we already have?
if (depth>=mindepth)
    return;
end;

stagethatoutputs=[];

% Search the current level, see if it provides the required output
if (~notthislevelplease)
    depth=depth+1;
    [stagepath stagename]=fileparts(aap.tasklist.main.module(currentstage).name);
    stagetag=aas_getstagetag(aap,currentstage);
    index=aap.tasklist.main.module(currentstage).index;
    
    if (isfield(aap.schema.tasksettings.(stagename)(index),'outputstreams'))
        outputstreams=aap.schema.tasksettings.(stagename)(index).outputstreams;
        for i=1:length(outputstreams.stream)
            if (strcmp(outputtype,outputstreams.stream{i}) || strcmp(outputtype,[stagetag '.' outputstreams.stream{i}]))
                stagethatoutputs=currentstage;
                mindepth=depth;
            end;
        end;
    end;
end;

% If not found, search backwards further
if (isempty(stagethatoutputs))
    dependenton=aap.internal.dependenton{currentstage};
    for i=1:length(dependenton)
        [aap stagethatoutputs mindepth]=aas_searchforoutput(aap,dependenton(i).stage,outputtype,false,depth,mindepth);
    end;
end;
