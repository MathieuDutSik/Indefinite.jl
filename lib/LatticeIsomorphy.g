FileISOM32:=Filename(DirectoriesPackagePrograms("indefinite"),"ISOM");
FileAUTO32:=Filename(DirectoriesPackagePrograms("indefinite"),"AUTO");
FileISOM64:=Filename(DirectoriesPackagePrograms("indefinite"),"ISOM");
FileAUTO64:=Filename(DirectoriesPackagePrograms("indefinite"),"AUTO");
FileISOMtoGAP:=Filename(DirectoriesPackagePrograms("indefinite"),"ISOMtoGAP");
FileAUTOMtoGAP:=Filename(DirectoriesPackagePrograms("indefinite"),"AUTOMtoGAP");


__ExtractInvariantZBasisShortVectorNoGroupGeneralLazy:=function(GramMat, GetSHVgroup, DoExpansion)
  local n, IsBetterProperty, IsWishedSet, ListVector, EWR, eSHV, posSHV, SHVgroup;
  n:=Length(GramMat);
  IsBetterProperty:=function(NewSet, OldSet)
    local rNew, rOld, detNew, detOld, eSelect, NewSetProj, OldSetProj;
    rNew:=PersoRankMat(NewSet);
    rOld:=PersoRankMat(OldSet);
    if rNew>rOld then
      return true;
    fi;
    if rNew=n then
      eSelect:=[1..n];
    else
      eSelect:=ColumnReduction(NewSet).Select;
    fi;
    NewSetProj:=List(NewSet, x->x{eSelect});
    OldSetProj:=List(OldSet, x->x{eSelect});
    detNew:=AbsInt(DeterminantMat(BaseIntMat(NewSetProj)));
    detOld:=AbsInt(DeterminantMat(BaseIntMat(OldSetProj)));
    if detNew<detOld then
      return true;
    fi;
    return false;
  end;
  IsWishedSet:=function(TheSet)
    local det;
    if PersoRankMat(TheSet)<n then
      return false;
    fi;
    det:=AbsInt(DeterminantMat(BaseIntMat(TheSet)));
    if det=1 then
      return true;
    fi;
    return false;
  end;
  ListVector:=[];
  SHVgroup:=[];
  posSHV:=0;
  while(true)
  do
    posSHV:=posSHV+1;
    if posSHV>Length(SHVgroup) then
      DoExpansion();
      SHVgroup:=GetSHVgroup();
    fi;
    eSHV:=SHVgroup[posSHV];
    EWR:=Concatenation(ListVector, eSHV);
    if IsBetterProperty(EWR, ListVector)=true then
      Append(ListVector, eSHV);
    fi;
    if IsWishedSet(ListVector)=true then
      if Length(ListVector)<>Length(Set(ListVector)) then
        Error("This is not correct");
      fi;
      return ListVector;
    fi;
  od;
end;


GetReducedDiagonal:=function(GramMat)
  local n, WorkGramMat, RedGramMat, TheDiag;
  n:=Length(GramMat);
  WorkGramMat:=RemoveFractionMatrix(GramMat);
  RedGramMat:=LLLReducedGramMat(WorkGramMat).remainder;
  TheDiag:=List([1..n], x->RedGramMat[x][x]);
  return TheDiag;
end;


__ExtractInvariantZBasisShortVectorNoGroup_V1:=function(GramMat)
  local n, WorkGramMat, LNorm, SHV, ListVector, ListNorms, eNorm, eSHV, iV, SHVgroup, TheDet, eSet, ListInv, eEnt, GetFctInvariant, TheColl, iNorm, DoExpansion, GetSHVgroup, TheDiag;
  if IsPositiveDefiniteSymmetricMatrix(GramMat)=false then
    Print("Error in __ExtractInvariantZBasisShortVectorNoGroup\n");
    Error("The matrix should be positive definite");
  fi;
  n:=Length(GramMat);
  WorkGramMat:=RemoveFractionMatrix(GramMat);
  TheDiag:=GetReducedDiagonal(GramMat);
#  Print("TheDiag=", TheDiag, "\n");
  LNorm:=Set(TheDiag);
  SHV:=ShortestVectors(WorkGramMat, Maximum(LNorm));
#  Print("|SHV.norms|=", Length(SHV.norms), " Coll=", Collected(SHV.norms), "\n");
  ListVector:=[];
  ListNorms:=Set(SHV.norms);
  TheDet:=Determinant(BaseIntMat(SHV.vectors));
  if AbsInt(TheDet)<>1 then
    Error("We have an inconsistency in determinant");
  fi;
  SHVgroup:=[];
  GetFctInvariant:=function(SHVgroup, eVect)
    return List(SHVgroup, y->Collected(List(y, x->x*WorkGramMat*eVect)));
  end;
  iNorm:=0;
  DoExpansion:=function()
    iNorm:=iNorm+1;
    eNorm:=ListNorms[iNorm];
#    Print("iNorm=", iNorm, " eNorm=", eNorm, "\n");
    eSHV:=[];
    for iV in [1..Length(SHV.norms)]
    do
      if SHV.norms[iV]=eNorm then
        Add(eSHV, SHV.vectors[iV]);
        Add(eSHV, -SHV.vectors[iV]);
      fi;
    od;
    if Length(SHVgroup)=1 then
      ListInv:=ListWithIdenticalEntries(Length(eSHV), eNorm);
    else
      ListInv:=List(eSHV, x->GetFctInvariant(SHVgroup, x));
    fi;
    TheColl:=Collected(ListInv);
    for eEnt in TheColl
    do
      eSet:=Filtered([1..Length(ListInv)], x->ListInv[x]=eEnt[1]);
      Add(SHVgroup, eSHV{eSet});
    od;
    Add(SHVgroup, eSHV);
  end;
  GetSHVgroup:=function()
    return SHVgroup;
  end;
  return __ExtractInvariantZBasisShortVectorNoGroupGeneralLazy(WorkGramMat, GetSHVgroup, DoExpansion);
end;


__ExtractInvariantZBasisShortVectorNoGroup:=__ExtractInvariantZBasisShortVectorNoGroup_V1;


DoesItContainStandardBasis:=function(SVR)
  local n, i, j, eSgn, TestVect;
  if Length(SVR)=0 then
    return false;
  fi;
  n:=Length(SVR[1]);
  for i in [1..n]
  do
    for j in [0..1]
    do
      eSgn:=2*j-1;
      TestVect:=ListWithIdenticalEntries(n, 0);
      TestVect[i]:=eSgn;
      if Position(SVR, TestVect)=fail then
        return false;
      fi;
    od;
  od;
  return true;
end;


AUTO_ISOM__PrintMatrix:=function(output, TheMat)
  local PrintLowMat, PrintTotalMat, n;
  n:=Length(TheMat);
  PrintLowMat:=function(TheMat)
    local i, j;
    for i in [1..n]
    do
      for j in [1..i]
      do
        AppendTo(output, " ", TheMat[i][j]);
      od;
      AppendTo(output, "\n");
    od;
  end;
  PrintTotalMat:=function(TheMat)
    local i, j;
    for i in [1..n]
    do
      for j in [1..n]
      do
        AppendTo(output, " ", TheMat[i][j]);
      od;
      AppendTo(output, "\n");
    od;
  end;
  if TheMat=TransposedMat(TheMat) then
    AppendTo(output, n, "x0\n");
    PrintLowMat(TheMat);
  elif TheMat=-TransposedMat(TheMat) then
    AppendTo(output, n, "x-1\n");
    PrintLowMat(TheMat);
  else
    AppendTo(output, n, "\n");
    PrintTotalMat(TheMat);
  fi;
end;


ArithmeticAutomorphismMatrixFamily_Souvignier:=function(TheOption, ListMat, SVR)
  local n, FileIn, FileOut, FileGap, FileError32, FileError64, output, eMat, eVect, i, TheNorm, ChainOption, REP, TheMatrixGRP, KeyTest, eGen, rListMat, TheMult, eRec, eLine, eEnt;
  n:=Length(ListMat[1]);
  if IsPositiveDefiniteSymmetricMatrix(ListMat[1])=false then
    Error("The first matrix should be positive definite");
  fi;
  TheMult:=1;
  for eMat in ListMat
  do
    for eLine in eMat
    do
      for eEnt in eLine
      do
        TheMult:=LcmInt(TheMult, DenominatorRat(eEnt));
      od;
    od;
  od;
  rListMat:=[];
  for eMat in ListMat
  do
    Add(rListMat, TheMult*eMat);
  od;
  FileIn:=Filename(POLYHEDRAL_tmpdir, "AUTO.in");
  FileOut:=Filename(POLYHEDRAL_tmpdir, "AUTO.out");
  FileGap:=Filename(POLYHEDRAL_tmpdir, "AUTO.gap");
  FileError32:=Filename(POLYHEDRAL_tmpdir, "AUTO.error32");
  FileError64:=Filename(POLYHEDRAL_tmpdir, "AUTO.error64");
  output:=OutputTextFile(FileIn, true);
  AppendTo(output, "#", Length(rListMat), "\n");
  for eMat in rListMat
  do
    AUTO_ISOM__PrintMatrix(output, eMat);
  od;
  ChainOption:="";
  KeyTest:=DoesItContainStandardBasis(SVR);
  if KeyTest=true then
    ChainOption:=" -V";
    AppendTo(output, "\n");
    AppendTo(output, Length(SVR), "\n");
    for eVect in SVR
    do
      for i in [1..n]
      do
        AppendTo(output, " ", eVect[i]);
      od;
      TheNorm:=eVect*rListMat[1]*eVect;
      AppendTo(output, " ", TheNorm, "\n");
    od;
  fi;
  CloseStream(output);
  # not completely satisfying here ...
  Exec(FileAUTO32, TheOption, ChainOption, " < ", FileIn, " > ", FileOut, "2>", FileError32);
  Exec(FileAUTOMtoGAP, " ", FileOut, " > ", FileGap);
  REP:=ReadAsFunction(FileGap)();
  TheMatrixGRP:=Group(REP.ListMat);
  SetSize(TheMatrixGRP, REP.order);
  RemoveFile(FileIn);
  RemoveFile(FileOut);
  RemoveFile(FileGap);
  RemoveFile(FileError32);
  RemoveFile(FileError64);
  # just a check to be sure
  for eGen in GeneratorsOfGroup(TheMatrixGRP)
  do
    for eMat in ListMat
    do
      if eGen*eMat*TransposedMat(eGen)<>eMat then
        Error("We have an error in ArithmeticAutomorphismMatrixFamily_Souvignier");
      fi;
    od;
    if Set(SVR*eGen)<>Set(SVR) then
      Print("The SVR is not stabilized. This is a feature, not a bug.\n");
      Print("But likely not what you want.\n");
    fi;
  od;
  return TheMatrixGRP;
end;


ArithmeticAutomorphismMatrixFamily_HackSouvignier:=function(TheOption, ListMat, SVR, AffBas)
  local ListMatMapped, SVRMapped, TheGRP, GRPret, eGen, eMat;
  if IsPositiveDefiniteSymmetricMatrix(ListMat[1])=false then
    Error("The first matrix should be positive definite");
  fi;
  ListMatMapped:=List(ListMat, x->AffBas*x*TransposedMat(AffBas));
  if Length(SVR)=0 then
    SVRMapped:=[];
  else
    SVRMapped:=SVR*Inverse(AffBas);
  fi;
  TheGRP:=ArithmeticAutomorphismMatrixFamily_Souvignier(TheOption, ListMatMapped, SVRMapped);
  GRPret:=Group(List(GeneratorsOfGroup(TheGRP), x->Inverse(AffBas)*x*AffBas));
  SetSize(GRPret, Size(TheGRP));
  for eGen in GeneratorsOfGroup(GRPret)
  do
    for eMat in ListMat
    do
      if eGen*eMat*TransposedMat(eGen)<>eMat then
        Error("Error in ArithmeticAutomorphismMatrixFamily_HackSouvignier");
      fi;
    od;
  od;
  return GRPret;
end;


LattIsom_CreateAffineBasis:=function(SVR)
  local ListCanStart, i, eVect, n;
  n:=Length(SVR[1]);
  ListCanStart:=[];
  for i in [1..n]
  do
    eVect:=ListWithIdenticalEntries(n,0);
    eVect[i]:=1;
    if Position(SVR, eVect)<>fail then
      Add(ListCanStart, eVect);
    fi;
  od;
  return ExtendToCompleteAffineBasis(SVR, ListCanStart);
end;


ArithmeticAutomorphismMatrixFamily_HackSouvignier_V2:=function(TheOption, ListMat, SVR)
  local TheP;
  if IsPositiveDefiniteSymmetricMatrix(ListMat[1])=false then
    Error("The first matrix should be positive definite");
  fi;
  TheP:=LattIsom_CreateAffineBasis(SVR);
  return ArithmeticAutomorphismMatrixFamily_HackSouvignier(TheOption, ListMat, SVR, TheP);
end;


ArithmeticAutomorphismGroup_InputBasis:=function(ListGramMat, InvariantBasis)
  return ArithmeticAutomorphismMatrixFamily_HackSouvignier_V2("", ListGramMat, InvariantBasis);
end;


ArithmeticAutomorphismGroup:=function(ListGramMat)
  local InvariantBasis;
  InvariantBasis:=__ExtractInvariantZBasisShortVectorNoGroup(ListGramMat[1]);
  return ArithmeticAutomorphismGroup_InputBasis(ListGramMat, InvariantBasis);
end;


ArithmeticEquivalenceMatrixFamily_Souvignier:=function(TheOption, ListMat1, SVR1, ListMat2, SVR2)
  local n, n1, n2, FileIn, FileOut, FileGap, FileError32, FileError64, output, ChainOption, eMat, eVect, i, TheNorm, REP, KeyTest, TheReply, iMat, TheClean, rListMat1, rListMat2, TheMult, eRec, eLine, eEnt;
  n1:=Length(ListMat1[1]);
  n2:=Length(ListMat2[1]);
  if n1<>n2 then
    Error("It is not allowed to have two input matrices of different size");
  fi;
  n:=n1;
  if Length(SVR1)=0 and Length(SVR2) > 0 then
    Error("Either you provide basis as input for both or neither");
  fi;
  if Length(SVR2)=0 and Length(SVR1) > 0 then
    Error("Either you provide basis as input for both or neither");
  fi;
  if IsPositiveDefiniteSymmetricMatrix(ListMat1[1])=false then
    Error("The first matrix of ListMat1 should be positive definite");
  fi;
  if IsPositiveDefiniteSymmetricMatrix(ListMat2[1])=false then
    Error("The first matrix of ListMat2 should be positive definite");
  fi;
  TheMult:=1;
  for eMat in ListMat1
  do
    for eLine in eMat
    do
      for eEnt in eLine
      do
        TheMult:=LcmInt(TheMult, DenominatorRat(eEnt));
      od;
    od;
  od;
  for eMat in ListMat2
  do
    for eLine in eMat
    do
      for eEnt in eLine
      do
        TheMult:=LcmInt(TheMult, DenominatorRat(eEnt));
      od;
    od;
  od;
  rListMat1:=[];
  for eMat in ListMat1
  do
    Add(rListMat1, TheMult*eMat);
  od;
  rListMat2:=[];
  for eMat in ListMat2
  do
    Add(rListMat2, TheMult*eMat);
  od;
  TheClean:=function()
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    RemoveFile(FileGap);
    RemoveFile(FileError32);
    RemoveFile(FileError64);
  end;
  FileIn:=Filename(POLYHEDRAL_tmpdir, "EQUIV.in");
  FileOut:=Filename(POLYHEDRAL_tmpdir, "EQUIV.out");
  FileGap:=Filename(POLYHEDRAL_tmpdir, "EQUIV.gap");
  FileError32:=Filename(POLYHEDRAL_tmpdir, "EQUIV.error32");
  FileError64:=Filename(POLYHEDRAL_tmpdir, "EQUIV.error64");
  output:=OutputTextFile(FileIn, true);
  AppendTo(output, "#", Length(ListMat1), "\n");
  for eMat in rListMat1
  do
    AUTO_ISOM__PrintMatrix(output, eMat);
  od;
  #
  KeyTest:=DoesItContainStandardBasis(SVR1) and DoesItContainStandardBasis(SVR2);
  ChainOption:="";
  if KeyTest=true then
    ChainOption:=" -V3";
    AppendTo(output, "\n");
    AppendTo(output, Length(SVR1), "\n");
    for eVect in SVR1
    do
      for i in [1..n]
      do
        AppendTo(output, " ", eVect[i]);
      od;
      TheNorm:=eVect*rListMat1[1]*eVect;
      AppendTo(output, " ", TheNorm, "\n");
    od;
  fi;
  for eMat in rListMat2
  do
    AUTO_ISOM__PrintMatrix(output, eMat);
  od;
  if KeyTest=true then
    AppendTo(output, "\n");
    AppendTo(output, Length(SVR2), "\n");
    for eVect in SVR2
    do
      for i in [1..n]
      do
        AppendTo(output, " ", eVect[i]);
      od;
      TheNorm:=eVect*rListMat2[1]*eVect;
      AppendTo(output, " ", TheNorm, "\n");
    od;
  fi;
  CloseStream(output);
  # not completely satisfying here ...
  Exec(FileISOM32, TheOption, ChainOption, " ", FileIn, " > ", FileOut, " 2>", FileError32);
  Exec(FileISOMtoGAP, " ", FileOut, " > ", FileGap);
  REP:=ReadAsFunction(FileGap)();
  if REP=false then
    TheClean();
    return false;
  fi;
  TheReply:=Inverse(REP);
  for iMat in [1..Length(ListMat1)]
  do
    if TheReply*ListMat1[iMat]*TransposedMat(TheReply)<>ListMat2[iMat] then
      Error("Error in ArithmeticEquivalenceMatrixFamily_Souvignier");
    fi;
  od;
  if Set(SVR1*Inverse(TheReply))<>Set(SVR2) then
    Print("The SVR1 / SVR2 are not mapped. This is a feature, not a bug.\n");
    Print("But likely not what you want.\n");
  fi;
  TheClean();
  return TheReply;
end;


ArithmeticEquivalenceMatrixFamily_HackSouvignier:=function(TheOption, ListMat1, SVR1, AffBas1, ListMat2, SVR2, AffBas2)
  local ListMat1mapped, SVR1mapped, ListMat2mapped, SVR2mapped, TheReplyMapped, TheReply, iMat;
  if IsPositiveDefiniteSymmetricMatrix(ListMat1[1])=false then
    Error("The first matrix of ListMat1 should be positive definite");
  fi;
  if IsPositiveDefiniteSymmetricMatrix(ListMat2[1])=false then
    Error("The first matrix of ListMat2 should be positive definite");
  fi;
  ListMat1mapped:=List(ListMat1, x->AffBas1*x*TransposedMat(AffBas1));
  if Length(SVR1)=0 then
    SVR1mapped:=[];
  else
    SVR1mapped:=SVR1*Inverse(AffBas1);
  fi;
  ListMat2mapped:=List(ListMat2, x->AffBas2*x*TransposedMat(AffBas2));
  if Length(SVR2)=0 then
    SVR2mapped:=[];
  else
    SVR2mapped:=SVR2*Inverse(AffBas2);
  fi;
  TheReplyMapped:=ArithmeticEquivalenceMatrixFamily_Souvignier(TheOption, ListMat1mapped, SVR1mapped, ListMat2mapped, SVR2mapped);
  if TheReplyMapped=false then
    return false;
  fi;
  TheReply:=Inverse(AffBas2)*TheReplyMapped*AffBas1;
  for iMat in [1..Length(ListMat1)]
  do
    if TheReply*ListMat1[iMat]*TransposedMat(TheReply)<>ListMat2[iMat] then
      Error("ArithmeticEquivalenceMatrixFamily_HackSouvignier: matrix inconsistency");
    fi;
  od;
  return TheReply;
end;


ArithmeticEquivalenceMatrixFamily_HackSouvignier_V2:=function(TheOption, ListMat1, SVR1, ListMat2, SVR2)
  local TheP1, TheP2;
  if IsPositiveDefiniteSymmetricMatrix(ListMat1[1])=false then
    Error("The first matrix of ListMat1 should be positive definite");
  fi;
  if IsPositiveDefiniteSymmetricMatrix(ListMat2[1])=false then
    Error("The first matrix of ListMat2 should be positive definite");
  fi;
  TheP1:=LattIsom_CreateAffineBasis(SVR1);
  TheP2:=LattIsom_CreateAffineBasis(SVR2);
  return ArithmeticEquivalenceMatrixFamily_HackSouvignier(TheOption, ListMat1, SVR1, TheP1, ListMat2, SVR2, TheP2);
end;


ArithmeticIsomorphism_Isom:=function(ListGramMat1, ListGramMat2)
  local InvariantBasis1, InvariantBasis2;
  InvariantBasis1:=__ExtractInvariantZBasisShortVectorNoGroup(ListGramMat1[1]);
  InvariantBasis2:=__ExtractInvariantZBasisShortVectorNoGroup(ListGramMat2[1]);
  return ArithmeticEquivalenceMatrixFamily_HackSouvignier_V2("", ListGramMat1, InvariantBasis1, ListGramMat2, InvariantBasis2);
end;


__CharacteristicGraphMatrixFamily:=function(ListMat, SVR)
  local nbV, DistMat, iVect, jVect, TheRec;
  nbV:=Length(SVR);
#  Print("nbV=", nbV, "\n");
#  Print("|ListMat|=", Length(ListMat), "\n");
  DistMat:=NullMat(nbV, nbV);
  for iVect in [1..nbV-1]
  do
    for jVect in [iVect+1..nbV]
    do
      TheRec:=List(ListMat, x->SVR[iVect]*x*SVR[jVect]);
      DistMat[iVect][jVect]:=TheRec;
      DistMat[jVect][iVect]:=TheRec;
    od;
  od;
  return DistMat;
end;


ArithmeticEquivalenceMatrixFamily_Nauty:=function(ListMat1, SVR1, ListMat2, SVR2)
  local DistMat1, DistMat2, eList, RetMat, iMat, RetMat2;
  if RankMat(SVR1)<>Length(ListMat1[1][1]) then
    Error("Wrong SVR1 argument");
  fi;
  if RankMat(SVR2)<>Length(ListMat2[1][1]) then
    Error("Wrong SVR2 argument");
  fi;
  DistMat1:=__CharacteristicGraphMatrixFamily(ListMat1, SVR1);
  DistMat2:=__CharacteristicGraphMatrixFamily(ListMat2, SVR2);
  eList:=IsIsomorphicEdgeColoredGraph(DistMat1, DistMat2);
  if eList=false then
    return false;
  fi;
  RetMat:=FindTransformation(SVR1, SVR2, PermList(eList));
  RetMat2:=Inverse(RetMat);
  for iMat in [1..Length(ListMat1)]
  do
    if RetMat2*ListMat1[iMat]*TransposedMat(RetMat2)<>ListMat2[iMat] then
      Error("Error in ArithmeticEquivalenceMatrixFamily_Nauty");
    fi;
  od;
  return RetMat2;
end;


ArithmeticIsomorphism_Nauty:=function(ListGramMat1, ListGramMat2)
  local InvariantBasis1, InvariantBasis2;
  InvariantBasis1:=__ExtractInvariantZBasisShortVectorNoGroup(ListGramMat1[1]);
  InvariantBasis2:=__ExtractInvariantZBasisShortVectorNoGroup(ListGramMat2[1]);
  return ArithmeticEquivalenceMatrixFamily_Nauty(ListGramMat1, InvariantBasis1, ListGramMat2, InvariantBasis2);
end;


ArithmeticIsomorphism:=function(ListGramMat1, ListGramMat2)
  if Length(ListGramMat1)<>Length(ListGramMat2) then
    return false;
  fi;
  if IsPositiveDefiniteSymmetricMatrix(ListGramMat1[1])=false then
    return false;
  fi;
  if IsPositiveDefiniteSymmetricMatrix(ListGramMat2[1])=false then
    return false;
  fi;
  return ArithmeticIsomorphism_Isom(ListGramMat1, ListGramMat2);
end;
