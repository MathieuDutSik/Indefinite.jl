RemoveFileIfExist:=function(FileName)
    local eFile;
    if IsString(FileName) then
        if IsExistingFile(FileName) then
            RemoveFile(FileName);
        fi;
    else
        if IsList(FileName) then
            for eFile in FileName
            do
                RemoveFileIfExist(eFile);
            od;
        fi;
    fi;
end;
