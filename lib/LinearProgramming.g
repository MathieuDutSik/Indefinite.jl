FileSamplingFacet:=Filename(DirectoriesPackagePrograms("indefinite"),"POLY_sampling_facets");
FileCddLp2:=Filename(DirectoriesPackagePrograms("indefinite"),"POLY_cdd_lp2");


#
# this thing is the end of a long design road.
# realizing that linear programming is quite complex.
# Basically, everything is outputed, everything is read
# and you have to make interpretations yourself.
LinearProgramming:=function(InequalitySet, ToBeMinimized)
  local FileIne, FileLps, FileErr, FileGap, FileDdl, FileLog, outputCdd, input, eLine, TheLP, TheDim, eVect, eSum, eEnt, nbIneq, TheCommand1, TheCommand2;
  FileIne:=Filename(POLYHEDRAL_tmpdir, "LP.listine");
  FileMin:=Filename(POLYHEDRAL_tmpdir, "LP.tominimize");
  FileRes:=Filename(POLYHEDRAL_tmpdir, "LP.result");
  RemoveFileIfExist(FileIne);
  RemoveFileIfExist(FileMin);
  RemoveFileIfExist(FileRes);
  WriteMatrixFile(FileIne, InequalitySet);
  WriteVectorFile(FileMin, ToBeMinimized);
  #
  TheCommand1:=Concatenation(FileCddLp2, " ", FileIne, " ", FileMin, " > ", FileRes);
  Exec(TheCommand1);
  #
  TheLP:=ReadAsFunction(FileRes)();
  RemoveFileIfExist(FileIne);
  RemoveFileIfExist(FileMin);
  RemoveFileIfExist(FileRes);
  return TheLP;
end;


GetInitialRaysGeneral:=function(FAC, command)
    local FileExt, FileOut;
    FileExt:=Concatenation(ThePath, "Sampling.ext");
    FileOut:=Concatenation(ThePath, "Sampling.out");
    Exec(FileDualDescriptionGroup, " rational ", command, " ", FileExt, " ", FileOut);
    ListInc:=ReadAsFunction(FileOutput)();
    if Length(ListInc)=0 then
        Error("Error in GetInitialRaysGeneral");
    fi;
    RemoveFile(FileExt);
    RemoveFile(FileOut);
    return ListIncd;
end;

GetInitialRays_LinProg:=function(EXT, nb)
    local command;
    command:=Concatenation("lp_cdd:iter_", String(nb));
    return GetInitialRaysGeneral(EXT, command);
end;


Norm_L1:=function(eVect)
  return Sum(List(eVect, AbsInt));
end;
