FileIsEmptyFile:=Filename(DirectoriesPackagePrograms("indefinite"),"IsEmptyFile");


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


IsEmptyFile:=function(FileName)
  local FileRep, TheRep;
  FileRep:=Filename(POLYHEDRAL_tmpdir, "FileInterpretation");
  Exec(FileIsEmptyFile, " ", FileName, " > ", FileRep);
  TheRep:=ReadAsFunction(FileRep)();
  RemoveFile(FileRep);
  return TheRep;
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


SaveDataToFileRecoverablePrevState:=function(FileName, Data)
  local File1, File2, File1touch, File2touch;
  File1:=Concatenation(FileName, "_1");
  File2:=Concatenation(FileName, "_2");
  File1touch:=Concatenation(File1, "_touch");
  File2touch:=Concatenation(File2, "_touch");
  SaveDataToFile(File2, Data);
  SaveDataToFile(File2touch, 0);
  RemoveFileIfExist(File1touch);
  RemoveFileIfExist(File1);
  Exec("cp ", File2, " ", File1);
  Exec("cp ", File2touch, " ", File1touch);
end;


SaveDataToFileRecoverablePrevStatePlusTest:=function(FileName, Data, test)
  if test then
    SaveDataToFileRecoverablePrevState(FileName, Data);
  fi;
end;


ReadAsFunctionRecoverablePrevState:=function(FileName)
  local File1, File2, File1touch, File2touch;
  File1:=Concatenation(FileName, "_1");
  File2:=Concatenation(FileName, "_2");
  File1touch:=Concatenation(File1, "_touch");
  File2touch:=Concatenation(File2, "_touch");
  if IsExistingFile(File1touch) then
    return ReadAsFunction(File1)();
  fi;
  if IsExistingFile(File2touch) then
    return ReadAsFunction(File2)();
  fi;
end;


IsExistingFileRecoverablePrevState:=function(FileName)
  local File1, File2, File1touch, File2touch;
  File1:=Concatenation(FileName, "_1");
  File2:=Concatenation(FileName, "_2");
  File1touch:=Concatenation(File1, "_touch");
  File2touch:=Concatenation(File2, "_touch");
  if IsExistingFile(File2touch) then
    return true;
  fi;
  if IsExistingFile(File1touch) then
    return true;
  fi;
  return false;
end;


SaveDebugInfo:=function(ePrefix, TheData)
    local n_index, FileSave;
    n_index:=0;
    while(true)
    do
        FileSave:=Concatenation(ePrefix, String(n_index));
	if IsExistingFile(FileSave)=false then
            SaveDataToFile(FileSave, TheData);
            return;
        fi;
        n_index:=n_index+1;
    od;
end;
