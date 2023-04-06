FileGLRS:=Filename(DirectoriesPackagePrograms("indefinite"),"glrs");
FileIsoReductionNG:=Filename(DirectoriesPackagePrograms("indefinite"),"POLY_IsomorphismReduction");
FileNudifyLRS_reductionNG:=Filename(DirectoriesPackagePrograms("indefinite"),"NudifyLRS.reductionNG");


__DualDescriptionLRS_Reduction:=function(EXT, GroupExt, ThePath)
  local eSub, EXT2, EXT3, FileExt, FileOut, FileFAC, FileGroup, FileSupport, FileOutput, FileError, output, DimEXT, test, EXTnew, ListInc;
#  Print("Entering polyhedral function LRS_Reduction |GRP|=", Order(GroupExt), "\n");
  FileExt:=Concatenation(ThePath, "LRS_Project.ext");
  FileOut:=Concatenation(ThePath, "LRS_Project.out");
  FileFAC:=Concatenation(ThePath, "LRS_Project.fac");
  FileGroup:=Concatenation(ThePath, "LRS_Project.group");
  FileSupport:=Concatenation(ThePath, "LRS_Project.supo");
  FileOutput:=Concatenation(ThePath, "LRS_Project.output");
  FileError:=Concatenation(ThePath, "LRS_Project.error");
  RemoveFileIfExist(FileExt);
  RemoveFileIfExist(FileOut);
  RemoveFileIfExist(FileFAC);
  RemoveFileIfExist(FileGroup);
  RemoveFileIfExist(FileSupport);
  RemoveFileIfExist(FileOutput);
  RemoveFileIfExist(FileError);
  #
  output:=OutputTextFile(FileExt, true);
  eSub:=__ProjectionFrame(EXT);
  EXT2:=List(EXT, x->x{eSub});
  EXT3:=List(EXT2, RemoveFraction);
  if TestConicness(EXT3) then
    EXTnew:=ShallowCopy(EXT3);
  else
    EXTnew:=List(EXT3, x->Concatenation([0], x));
  fi;
  DimEXT:=Length(EXTnew[1]);
  #
  AppendTo(output, "V-representation\n");
  AppendTo(output, "begin\n");
  AppendTo(output, Length(EXTnew), " ", DimEXT, " integer\n");
  WriteMatrix(output, EXTnew);
  AppendTo(output, "end\n");
  CloseStream(output);
  #
  Exec(FileGLRS, " ", FileExt, " > ", FileOut);
  Exec(FileNudifyLRS_reductionNG, " ", FileFAC, " < ", FileOut);
  #
  WriteMatrixFile(FileSupport, EXTnew);
  #
  SYMPOL_PrintGroup(FileGroup, Length(EXTnew), GroupExt);
  #
  Exec(FileIsoReductionNG, " ", FileSupport, " ", FileFAC, " ", FileGroup, " ", FileOutput, "2>", FileError);
  ListInc:=ReadAsFunction(FileOutput)();
  if Length(ListInc)=0 then
    Error("Error in DualDescriptionLRS_Reduction");
  fi;
  RemoveFile(FileExt);
  RemoveFile(FileOut);
  RemoveFile(FileFAC);
  RemoveFile(FileGroup);
  RemoveFile(FileSupport);
  RemoveFile(FileOutput);
  RemoveFile(FileError);
  return ListInc;
end;
