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


IsExistingFilePlusTouch:=function(FileName)
  local FileTouch;
  if IsExistingFile(FileName)=false then
    return false;
  fi;
  FileTouch:=Concatenation(FileName, "_touch");
  return IsExistingFile(FileTouch);
end;


RemoveFileIfExistPlusTouch:=function(FileName)
  local FileTouch;
  RemoveFileIfExist(FileName);
  FileTouch:=Concatenation(FileName, "_touch");
  RemoveFileIfExist(FileTouch);
end;


RemoveFileIfExistPlusTouchPlusTest:=function(FileName, test)
  if test then
    RemoveFileIfExistPlusTouch(FileName);
  fi;
end;


IsExistingFilePlusTouchPlusTest:=function(FileName, test)
  local FileTouch;
  if test=false then
    return false;
  else
    FileTouch:=Concatenation(FileName, "_touch");
    if IsExistingFile(FileName)=false then
      return false;
    fi;
    return IsExistingFile(FileTouch);
  fi;
end;
