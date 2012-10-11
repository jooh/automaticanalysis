%% Recursively copy across aap parameters from structure
function destaap=aas_copyparameters(srcaap,destaap,nme)

fn=fieldnames(srcaap);
for fnind=1:length(fn)
    if (isstruct(srcaap.(fn{fnind})))
        % handle structure array parameters
        destlen = length(destaap.(fn{fnind}));
        srclen = length(srcaap.(fn{fnind}));
        % if default array already had a length > srclen something is
        % probably wrong
        assert(destlen<=srclen,'mismatched struct arrays')
        % upcast destlen to srclen
        if destlen<srclen
            destaap.(fn{fnind}) = destaap.(fn{fnind})(ones(1,srclen));
        end
        % iterate over struct array entries
        for n = 1:length(srcaap.(fn{fnind}))
            destaap.(fn{fnind})(n) = aas_copyparameters(...
                srcaap.(fn{fnind})(n),destaap.(fn{fnind})(n),...
                [nme '.' fn{fnind}]);
        end
    else
        if (~isfield(destaap,fn{fnind}))
            aas_log(srcaap,true,sprintf('Error when copying extra parameters, field %s is present in %s of extraparameters.aap  but not in normal aap structure',fn{fnind},nme));
        end;
        destaap.(fn{fnind})=srcaap.(fn{fnind});
    end;
end;
end
