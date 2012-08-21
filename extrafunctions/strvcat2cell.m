function strCell = strvcat2cell(strVcat)
    if size(strVcat,1)==1
        strCell = {strVcat};
    else
        strCell = {};
        for o = 1:size(strVcat,1)
            strCell = [strCell deblank(strVcat(o,:))];
        end
    end
end