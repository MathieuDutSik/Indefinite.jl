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


SaveDataToFile:=function(FileName, OBJ)
  local output;
  Exec("rm -f ", FileName,"\n");
  output:=OutputTextFile(FileName, true);;
  AppendTo(output, "return ", OBJ, ";\n");
  CloseStream(output);
end;


SaveDataToFilePlusTouch:=function(FileName, OBJ)
  local FileTouch;
  FileTouch:=Concatenation(FileName, "_touch");
  RemoveFileIfExist(FileTouch);
  SaveDataToFile(FileName, OBJ);
  SaveDataToFile(FileTouch, 0);
end;


IsExistingFilePlusTouch:=function(FileName)
  local FileTouch;
  if IsExistingFile(FileName)=false then
    return false;
  fi;
  FileTouch:=Concatenation(FileName, "_touch");
  return IsExistingFile(FileTouch);
end;


SaveDataToFilePlusTouchPlusTest:=function(FileName, OBJ, test)
  if test then
    SaveDataToFilePlusTouch(FileName, OBJ);
  fi;
end;


RemoveDirectoryPlusTest:=function(FileDirectory, test)
  if test then
    Exec("rm -rf ", FileDirectory);
  fi;
end;


CreateDirectoryPlusTest:=function(FileDirectory, test)
  if test then
    Exec("mkdir -p ", FileDirectory);
  fi;
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


ComputeAndSavePlusTouch:=function(FileName, FCT)
  local TheData;
  if IsExistingFilePlusTouch(FileName) then
    return ReadAsFunction(FileName)();
  else
    TheData:=FCT(1);
    SaveDataToFilePlusTouch(FileName, TheData);
    return TheData;
  fi;
end;


ComputeAndSavePlusTouchPlusTest:=function(FileName, FCT, test)
  if test=false then
    return FCT(1);
  else
    return ComputeAndSavePlusTouch(FileName, FCT);
  fi;
end;
