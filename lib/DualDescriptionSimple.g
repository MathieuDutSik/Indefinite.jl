FileDualDescriptionGroup:=Filename(DirectoriesPackagePrograms("indefinite"),"POLY_dual_description_group");


__DualDescriptionLRS_Reduction:=function(EXT, GroupExt, ThePath)
  local eSub, EXT2, EXT3, FileExt, FileOut, FileFAC, FileGroup, FileSupport, FileOutput, FileError, output, DimEXT, test, EXTnew, ListInc;
#  Print("Entering polyhedral function LRS_Reduction |GRP|=", Order(GroupExt), "\n");
  FileExt:=Concatenation(ThePath, "LRS_Project.ext");
  FileGrp:=Concatenation(ThePath, "LRS_Project.grp");
  FileOut:=Concatenation(ThePath, "LRS_Project.out");
  RemoveFileIfExist(FileExt);
  RemoveFileIfExist(FileGrp);
  RemoveFileIfExist(FileOut);
  #
  output:=OutputTextFile(FileExt, true);
  eSub:=__ProjectionFrame(EXT);
  EXT2:=List(EXT, x->x{eSub});
  WriteMatrixFile(FileSupport, EXTnew);
  SYMPOL_PrintGroup(FileGroup, Length(EXT2), GroupExt);
  #
  Exec(FileDualDescriptionGroup, " rational lrs_ring ", FileExt, " ", FileGrp, " ", FileOut);
  ListInc:=ReadAsFunction(FileOutput)();
  if Length(ListInc)=0 then
    Error("Error in DualDescriptionLRS_Reduction");
  fi;
  RemoveFile(FileExt);
  RemoveFile(FileGrp);
  RemoveFile(FileOut);
  return ListInc;
end;
