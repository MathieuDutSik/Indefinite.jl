Filelpcddcleaner:=Filename(DirectoriesPackagePrograms("indefinite"),"lpcddcleaner");
FileTestlp2:=Filename(DirectoriesPackagePrograms("indefinite"),"testlp2_gmp");


#
# this thing is the end of a long design road.
# realizing that linear programming is quite complex.
# Basically, everything is outputed, everything is read
# and you have to make interpretations yourself.
LinearProgramming_Rational:=function(InequalitySet, ToBeMinimized)
  local FileIne, FileLps, FileErr, FileGap, FileDdl, FileLog, outputCdd, input, eLine, TheLP, TheDim, eVect, eSum, eEnt, nbIneq, TheCommand1, TheCommand2;
  FileIne:=Filename(POLYHEDRAL_tmpdir, "LP.ine");
  FileLps:=Filename(POLYHEDRAL_tmpdir, "LP.lps");
  FileErr:=Filename(POLYHEDRAL_tmpdir, "LP.error");
  FileGap:=Filename(POLYHEDRAL_tmpdir, "LP.gap");
  FileDdl:=Filename(POLYHEDRAL_tmpdir, "LP.ddl");
  FileLog:=Filename(POLYHEDRAL_tmpdir, "LP.log");
#  Print("FileIne=", FileIne, "\n");
  RemoveFileIfExist(FileIne);
  RemoveFileIfExist(FileLps);
  RemoveFileIfExist(FileErr);
  RemoveFileIfExist(FileGap);
  RemoveFileIfExist(FileDdl);
  RemoveFileIfExist(FileLog);
  TheDim:=Length(InequalitySet[1]);
  nbIneq:=Length(InequalitySet);
#  Print("nbIneq=", nbIneq, "  TheDim=", TheDim, "\n");
  for eVect in InequalitySet
  do
    if Length(eVect)<>TheDim then
      Print("Incoherence in dimensions of InequalitySet, please");
      Error("Please correct");
    fi;
  od;
  if TheDim<>Length(ToBeMinimized) then
    Error("Incoherence in dimensions, please be careful");
  fi;
  outputCdd:=OutputTextFile(FileIne, true);;
  AppendTo(outputCdd, "H-representation\n");
  AppendTo(outputCdd, "begin\n");
  AppendTo(outputCdd, " ", Length(InequalitySet), " ", Length(ToBeMinimized), " integer\n");
  WriteMatrix(outputCdd, InequalitySet);
  AppendTo(outputCdd, "end\n");
  AppendTo(outputCdd, "minimize\n");
  WriteVector(outputCdd, ToBeMinimized);
  CloseStream(outputCdd);
  #
  TheCommand1:=Concatenation(FileTestlp2, " ", FileIne, " 2> ", FileErr, " > ", FileLog);
#  Print("TheCommand1=", TheCommand1, "\n");
  Exec(TheCommand1);
  #
  TheCommand2:=Concatenation(Filelpcddcleaner, " < ", FileLog, " > ", FileGap);
#  Print("TheCommand2=", TheCommand2, "\n");
  Exec(TheCommand2);
  #
  TheLP:=ReadAsFunction(FileGap)();
  if TheLP=rec() then
      Error("Debug from that point. TheLP=rec() is the error");
  fi;
  TheLP.method:="cdd";
  if IsBound(TheLP.dual_direction) then
#    Print("TheLP.dual_direction=", TheLP.dual_direction, "\n");
    eSum:=ListWithIdenticalEntries(TheDim,0);
    for eEnt in TheLP.dual_direction
    do
#      Print("eIneq=", InequalitySet[eEnt[1]], "\n");
      eSum:=eSum+InequalitySet[eEnt[1]]*AbsInt(eEnt[2]);
    od;
    if eSum[1] >= 0 then
      Print("Apparently something is not understood for\n");
      Print("cdd and linear programming (unfeasibilities) 1\n");
      Error("Please correct");
    fi;
    if eSum{[2..TheDim]}<>ListWithIdenticalEntries(TheDim-1,0) then
      Print("Apparently something is not understood for\n");
      Print("cdd and linear programming (unfeasibilities) 2\n");
      Error("Please correct");
    fi;
  fi;
  if IsBound(TheLP.primal_solution) then
    eVect:=ListWithIdenticalEntries(TheDim-1,0);
    for eEnt in TheLP.primal_solution
    do
      eVect[eEnt[1]]:=eEnt[2];
    od;
    TheLP.eVect:=eVect;
    TheLP.TheVert:=Concatenation([1], eVect);
  fi;
#  Print("Copy the files here\n");
#  Print(NullMat(5));
  RemoveFileIfExist(FileIne);
  RemoveFileIfExist(FileLps);
  RemoveFileIfExist(FileErr);
  RemoveFileIfExist(FileGap);
  RemoveFileIfExist(FileDdl);
  RemoveFileIfExist(FileLog);
  return TheLP;
end;


LinearProgramming:=function(InequalitySet, ToBeMinimized)
  local Nval;
  if IsMatrixRational(InequalitySet) and IsVectorRational(ToBeMinimized) then
    return LinearProgramming_Rational(InequalitySet, ToBeMinimized);
  fi;
  Error("You have to build your own arithmetic or use LinearProgrammingGeneralCode");
end;


GetPolytopizationInfo:=function(FAC)
  local eSet, FACred, dimRed, ListIneq, ToBeMinimized, eFac, eIneq, TheLP, eVect, eEnt, eMatBig, eBasis, eMatTrans, NewListIneq, eProd, fIneq;
  eSet:=ColumnReduction(FAC).Select;
  FACred:=List(FAC, x->x{eSet});
  dimRed:=Length(eSet);
  ListIneq:=[];
  ToBeMinimized:=ListWithIdenticalEntries(1+dimRed,0);
  for eFac in FACred
  do
    eIneq:=Concatenation([-1], eFac);
    Add(ListIneq, eIneq);
    ToBeMinimized:=ToBeMinimized+eIneq;
  od;
  TheLP:=LinearProgramming(ListIneq, ToBeMinimized);
  if IsBound(TheLP.primal_solution)=false then
    Error("Clear inconsistency, maybe the cone is empty");
  fi;
  eVect:=ListWithIdenticalEntries(dimRed,0);
  for eEnt in TheLP.primal_solution
  do
    eVect[eEnt[1]]:=eEnt[2];
  od;
  eMatBig:=Concatenation([eVect], IdentityMat(dimRed));
  eBasis:=RowReduction(eMatBig).EXT;
  eMatTrans:=TransposedMat(eBasis);
  return rec(eSet:=eSet, eMatTrans:=eMatTrans);
end;


DoPolytopization:=function(eRecPoly, FAC)
  local NewListIneq, FACred, eIneq, eProd, fIneq;
  NewListIneq:=[];
  FACred:=List(FAC, x->x{eRecPoly.eSet});
  for eIneq in FACred
  do
    eProd:=eIneq*eRecPoly.eMatTrans;
    fIneq:=eProd/eProd[1];
    Add(NewListIneq, fIneq);
  od;
  return NewListIneq;
end;


PolytopizationGeneralCone:=function(FAC)
  local eRecPoly;
  eRecPoly:=GetPolytopizationInfo(FAC);
  return DoPolytopization(eRecPoly, FAC);
end;


LP_GetSet:=function(FAC, eMinimize)
  local TheLP, eVect, eEnt, eVectExt, eSet, nbVert, TheDim;
  TheLP:=LinearProgramming(FAC, eMinimize);
  nbVert:=Length(FAC);
  TheDim:=Length(FAC[1]);
  if IsBound(TheLP.primal_solution) then
    eVect:=ListWithIdenticalEntries(TheDim,0);
    for eEnt in TheLP.primal_solution
    do
      eVect[eEnt[1]]:=eEnt[2];
    od;
    eVectExt:=Concatenation([1], eVect);
    eSet:=Filtered([1..nbVert], x->FAC[x]*eVectExt=0);
    if RankMat(FAC{eSet})=TheDim-1 then
      return eSet;
    fi;
    return fail;
  fi;
  return fail;
end;


LP_GetExpressionForLP:=function(EXT)
  local OneVert, ListDiff, eEXT, TheBasis, TheDim, NewListCoord, eVectSol, eIso, eVect;
  OneVert:=EXT[1];
  ListDiff:=[];
  for eEXT in EXT
  do
    Add(ListDiff, eEXT-OneVert);
  od;
  TheBasis:=RowReduction(ListDiff).EXT;
  TheDim:=Length(TheBasis);
  NewListCoord:=[];
  for eVect in ListDiff
  do
    eVectSol:=Concatenation([1], SolutionMat(TheBasis, eVect));
    Add(NewListCoord, eVectSol);
  od;
  eIso:=Isobarycenter(NewListCoord);
  eIso[1]:=0;
  return List(NewListCoord, x->x-eIso);
end;


GetInitialRaysEXT:=function(EXT, nb)
  local nbVert, EXT_lp, GetSet, ListSet, i, TheDim;
  nbVert:=Length(EXT);
  if Set(List(EXT, x->x[1]))<>[1] then
    Error("The first coordinate needs to be 1 for all vectors");
  fi;
  EXT_lp:=LP_GetExpressionForLP(EXT);
  TheDim:=Length(EXT_lp[1])-1;
  GetSet:=function()
    local siz, eMinimize, i, test;
    siz:=10;
    while(true)
    do
      eMinimize:=[1];
      for i in [1..TheDim]
      do
        Add(eMinimize, Random([-siz..siz]));
      od;
      test:=LP_GetSet(EXT_lp, eMinimize);
      if test<>fail then
        return test;
      fi;
      siz:=siz+1;
    od;
  end;
  ListSet:=[];
  for i in [1..nb]
  do
#    Print("GetInitialRaysEXT i=", i, " / ", nb, "\n");
    Add(ListSet, GetSet());
  od;
  return ListSet;
end;


GetInitialRaysGeneral:=function(FAC, nb)
  local NewListIneq;
  NewListIneq:=PolytopizationGeneralCone(FAC);
  return GetInitialRaysEXT(NewListIneq, nb);
end;


GetInitialRays_LinProg:=function(EXT, nb)
  local IsEXT, eEXT;
  Print("Beginning of GetInitialRays_LinProg\n");
  IsEXT:=true;
  for eEXT in EXT
  do
    if eEXT[1]<>1 then
      IsEXT:=false;
    fi;
  od;
  Print("IsEXT=", IsEXT, "\n");
  if IsEXT=true then
    return GetInitialRaysEXT(EXT, nb);
  fi;
  return GetInitialRaysGeneral(EXT, nb);
end;


Norm_L1:=function(eVect)
  return Sum(List(eVect, AbsInt));
end;
