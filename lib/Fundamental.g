FileGetDate:=Filename(DirectoriesPackagePrograms("indefinite"),"GetDate");


MatrixToVector:=function(eMat)
    local TheRet, dim, i;
    TheRet:=[];
    dim:=Length(eMat);
    for i in [1..dim]
    do
        Append(TheRet, eMat[i]);
    od;
    return TheRet;
end;


VectorToMatrix:=function(eVect)
    local dimSqr, dim, eMat, pos, i, j;
    dimSqr:=Length(eVect);
    dim:=Sqrt(dimSqr);
    eMat:=NullMat(dim,dim);
    pos:=0;
    for i in [1..dim]
    do
        for j in [1..dim]
        do
            pos:=pos+1;
            eMat[i][j]:=eVect[pos];
        od;
    od;
    return eMat;
end;


IntegralSpaceSaturation:=function(TheSpace)
    local dim, rnk, eOrth, eOrthB, TheSpaceInt;
    if Length(TheSpace)=0 then
        return [];
    fi;
    dim:=Length(TheSpace[1]);
    rnk:=RankMat(TheSpace);
    if rnk = dim then
        return IdentityMat(dim);
    fi;
    if rnk = 0 then
        return [];
    fi;
    eOrth:=NullspaceMat(TransposedMat(TheSpace));
    eOrthB:=List(eOrth, RemoveFraction);
    TheSpaceInt:=NullspaceIntMat(TransposedMat(eOrthB));
    return TheSpaceInt;
end;


#
#
# larger integer smaller than eFrac
LowerInteger:=function(eFrac)
  local a, b, r;
  if IsInt(eFrac)=true then
    return eFrac;
  fi;
  a:=NumeratorRat(eFrac);
  b:=DenominatorRat(eFrac);
  r:=a mod b;
  return (a-r)/b;
end;


# this procedure Build the Set:  Seto1 x Seto2 x .... x SetoN
BuildSetMultiple:=function(ListSeto)
  local DO, eSet, C, eC, eVal, U;
  DO:=[[]];
  for eSet in ListSeto
  do
    C:=ShallowCopy(DO);
    DO:=ShallowCopy([]);
    for eC in C
    do
      for eVal in eSet
      do
	U:=ShallowCopy(eC);
        Add(U, eVal);
        Add(DO, U);
      od;
    od;
  od;
  return DO;
end;


Isobarycenter:=function(EXT)
  return Sum(EXT)/Length(EXT);
end;


IsVectorRational:=function(eVect)
  local eX;
  for eX in eVect
  do
    if IsRat(eX)=false then
      return false;
    fi;
  od;
  return true;
end;


IsMatrixRational:=function(eMat)
  local eVect;
  for eVect in eMat
  do
    if IsVectorRational(eVect)=false then
      return false;
    fi;
  od;
  return true;
end;


TestConicness:=function(FAC)
  local eVal;
  for eVal in FAC
  do
    if eVal[1]<>0 then
      return false;
    fi;
  od;
  return true;
end;


FacetOfInfinity:=function(n)
  local VZ;
  VZ:=ListWithIdenticalEntries(n, 0);
  VZ[1]:=1;
  return VZ;
end;


IsIntegralVector:=function(eVect)
  local eVal;
  for eVal in eVect
  do
    if IsInt(eVal)=false then
      return false;
    fi;
  od;
  return true;
end;


IsIntegralMat:=function(eMat)
  local eLine, eVal;
  for eLine in eMat
  do
    for eVal in eLine
    do
      if IsInt(eVal)=false then
        return false;
      fi;
    od;
  od;
  return true;
end;


GetDate:=function()
  local FileD1, FileD2, reply;
  FileD1:=Filename(POLYHEDRAL_tmpdir,"DateFile1");
  FileD2:=Filename(POLYHEDRAL_tmpdir,"DateFile2");
  Exec("date +%s > ", FileD1);
  Exec(FileGetDate, " ", FileD1, " > ", FileD2);
  reply:=ReadAsFunction(FileD2)();
  RemoveFile(FileD1);
  RemoveFile(FileD2);
  return reply;
end;


RandomSubset:=function(eSet, k)
  local i, sSet, V, h;
  sSet:=[];
  V:=ListWithIdenticalEntries(Length(eSet), 1);
  for i in [1..k]
  do
    while(true)
    do
      h:=Random([1..Length(eSet)]);
      if V[h]=1 then
        V[h]:=0;
        Add(sSet, eSet[h]);
        break;
      fi;
    od;
  od;
  return Set(sSet);
end;


WriteVector:=function(outputarg, eLine)
  local eVal;
  for eVal in eLine
  do
    WriteAll(outputarg, Concatenation(" ", String(eVal)));
  od;
  WriteAll(outputarg, "\n");
end;


WriteMatrix:=function(outputarg, eMat)
  local eEXT;
  for eEXT in eMat
  do
    WriteVector(outputarg, eEXT);
  od;
end;


WriteMatrixFile:=function(eFile, EXT)
    local output, eEXT;
    output:=OutputTextFile(eFile, true);
    AppendTo(output, Length(EXT), " ", Length(EXT[1]), "\n");
    for eEXT in EXT
    do
        WriteVector(output, eEXT);
    od;
    CloseStream(output);
end;


CPP_WriteMatrix:=function(output, eMat)
  local nbRow, nbCol, iRow, iCol;
  nbRow:=Length(eMat);
  nbCol:=Length(eMat[1]);
  AppendTo(output, nbRow, " ", nbCol, "\n");
  for iRow in [1..nbRow]
  do
    for iCol in [1..nbCol]
    do
      AppendTo(output, " ", eMat[iRow][iCol]);
    od;
    AppendTo(output, "\n");
  od;
end;


RowReduction:=function(arg)
  local MatrixWork, rankMatrixWork, PreSelect, Irr, rank, pos, TE, SelectSet;
  if IsBound(arg[4]) then
    Error("Too many argument to this function");
  elif IsBound(arg[3]) then
    MatrixWork:=arg[1];
    rankMatrixWork:=arg[2];
    PreSelect:=arg[3];
  elif IsBound(arg[2]) then
    MatrixWork:=arg[1];
    rankMatrixWork:=arg[2];
    PreSelect:=[];
  elif IsBound(arg[1]) then
    MatrixWork:=arg[1];
    rankMatrixWork:=RankMat(MatrixWork);
    PreSelect:=[];
  else
    Error("Please provide at least one argument");
  fi;
  if Length(PreSelect)>0 then
    if RankMat(MatrixWork{PreSelect}) < Length(PreSelect) then
      Error("your preselection will not work to the task");
    else
      SelectSet:=ShallowCopy(PreSelect);
      Irr:=MatrixWork{SelectSet};
      rank:=Length(PreSelect);
    fi;
  else
    Irr:=[];
    SelectSet:=[];
    rank:=0;
  fi;
  pos:=1;
  SelectSet:=[];
  while (rank<rankMatrixWork)
  do
    if Position(SelectSet, pos)=fail then
      TE:=ShallowCopy(Irr);
      Add(TE, MatrixWork[pos]);
      if (RankMat(TE) > rank) then
        Add(Irr, MatrixWork[pos]);
        Add(SelectSet, pos);
        rank:=rank+1;
      fi;
    fi;
    pos:=pos+1;
  od;
  return rec(EXT:=Irr, Select:=SelectSet);
end;


ColumnReduction:=function(arg)
  local MatrixWork, rankMatrixWork, PreSelect, rep;
  if IsBound(arg[4]) then
    Error("Too many argument to this function");
  elif IsBound(arg[3]) then
    MatrixWork:=arg[1];
    rankMatrixWork:=arg[2];
    PreSelect:=arg[3];
  elif IsBound(arg[2]) then
    MatrixWork:=arg[1];
    rankMatrixWork:=arg[2];
    PreSelect:=[];
  elif IsBound(arg[1]) then
    MatrixWork:=arg[1];
    rankMatrixWork:=RankMat(MatrixWork);
    PreSelect:=[];
  else
    Error("Please provide at least one argument");
  fi;
  rep:=RowReduction(TransposedMat(MatrixWork), rankMatrixWork, PreSelect);
  return rec(EXT:=TransposedMat(rep.EXT), Select:=rep.Select);
end;


GcdVector:=function(TheVector)
  local eVal, TheVectorRed, LSred, LS1, ListCoef;
  if Length(TheVector)=1 then
    eVal:=TheVector[1];
    if eVal > 0 then
      return rec(TheGcd:=eVal, ListCoef:=[1]);
    fi;
    if eVal < 0 then
      return rec(TheGcd:=-eVal, ListCoef:=[-1]);
    fi;
    return rec(TheGcd:=0, ListCoef:=[1]);
  fi;
  TheVectorRed:=TheVector{[2..Length(TheVector)]};
  LSred:=GcdVector(TheVectorRed);
  LS1:=Gcdex(TheVector[1], LSred.TheGcd);
  ListCoef:=Concatenation([LS1.coeff1], LS1.coeff2*LSred.ListCoef);
  if ListCoef * TheVector <> LS1.gcd then
      Error("The family of coefficient does not get us the gcd");
  fi;
  return rec(TheGcd:=LS1.gcd, ListCoef:=ListCoef);
end;


GetZbasis:=function(ListElements)
  local TheDim, ListEqua, TheBasis, InvMatrix, eSet, GetOneBasis, ComputeSpeedingElements, FuncInsert, eElt, eEltRed, eSol, eLine, TheMult, DoCheck;
  if Length(ListElements)=0 then
    Print("|ListElements|=", Length(ListElements), "\n");
    Error("Problem in GetZbasis, we need at least one element");
  fi;
  TheDim:=Length(ListElements[1]);
  if RankMat(ListElements)=0 then
    return [];
  fi;
  ListEqua:=IdentityMat(TheDim);
  TheMult:=RemoveFractionMatrixPlusCoef(ListElements).TheMult;
  TheBasis:=[];
  InvMatrix:="irrelevant";
  eSet:=[];
  GetOneBasis:=function(eSol)
    local DimLoc, TheRedMat, eVect, iCol, n2, ListIdx, AbsList, TheMin, pos, ThePivot, r, q, NSP;
    DimLoc:=Length(TheBasis);
    TheRedMat:=Concatenation(IdentityMat(DimLoc), [eSol]);
    NSP:=NullspaceIntMat(RemoveFractionMatrix(TheRedMat));
    if Length(NSP)<>1 then
      Error("problem in computation of relation");
    fi;
    eVect:=ShallowCopy(NSP[1]);
    n2:=Length(eVect);
    while(true)
    do
      ListIdx:=Filtered([1..n2], x->eVect[x]<>0);
      if Length(ListIdx)=1 then
        return TheRedMat{Difference([1..n2], ListIdx)};
      fi;
      AbsList:=List(eVect{ListIdx}, AbsInt);
      TheMin:=Minimum(AbsList);
      pos:=Position(AbsList, TheMin);
      ThePivot:=ListIdx[pos];
      for iCol in Difference([1..n2], [ThePivot])
      do
        r:=eVect[iCol] mod AbsInt(eVect[ThePivot]);
        q:=(eVect[iCol] - r)/eVect[ThePivot];
        TheRedMat[ThePivot]:=TheRedMat[ThePivot]+q*TheRedMat[iCol];
        eVect[iCol]:=r;
        if First(eVect*TheRedMat, x->x<>0)<>fail then
          Error("inconsistency");
        fi;
      od;
    od;
  end;
  ComputeSpeedingElements:=function()
    ListEqua:=NullspaceMat(TransposedMat(TheBasis));
    eSet:=ColumnReduction(TheBasis).Select;
    InvMatrix:=Inverse(List(TheBasis, x->x{eSet}));
  end;
  FuncInsert:=function(eElt)
    local ListScal, TheSum, eEltRed, eSol, NewBasis, TheBasisNew;
    ListScal:=List(ListEqua, x->x*eElt);
    TheSum:=Sum(List(ListScal, x->x*x));
    if TheSum>0 then
      Add(TheBasis, eElt);
      TheBasis:=LLLReducedBasis(TheBasis).basis;
      ComputeSpeedingElements();
    else
      if Length(TheBasis)=0 then
        return;
      fi;
      eEltRed:=eElt{eSet};
      eSol:=eEltRed*InvMatrix;
      if IsIntegralVector(eSol) then
        return;
      fi;
      NewBasis:=GetOneBasis(eSol);
      TheBasisNew:=NewBasis*TheBasis;
      TheBasis:=LLLReducedBasis(TheBasisNew).basis;
      ComputeSpeedingElements();
    fi;
  end;
  for eElt in ListElements
  do
    FuncInsert(eElt);
  od;
  for eElt in ListElements
  do
    eEltRed:=eElt{eSet};
    eSol:=eEltRed*InvMatrix;
    if IsIntegralVector(eSol)=false then
      Error("Inconsistency in basis computation");
    fi;
  od;
  DoCheck:=true;
  if DoCheck then
    for eLine in TheBasis
    do
      eSol:=SolutionIntMat(ListElements*TheMult, eLine*TheMult);
      if eSol=fail then
        Error("Error in GetZbasis 1");
      fi;
    od;
    for eLine in ListElements
    do
      eSol:=SolutionIntMat(TheBasis*TheMult, eLine*TheMult);
      if eSol=fail then
        Error("Error in GetZbasis 2");
      fi;
    od;
  fi;
  return TheBasis;
end;


GetZbasisColumn:=function(Amat)
    return TransposedMat(GetZbasis(TransposedMat(Amat)));
end;


FindTransformation:=function(ListVert1, ListVert2, eGen)
  local eMat1, eMat2, iVert, jVert, rnk, Irr1, Irr2, rank, pos, TE, TheMat, test1, test2;
  eMat1:=[];
  eMat2:=[];
  if Length(ListVert1)<>Length(ListVert2) then
    Error("Error in FindTransformation: vector length is not the same");
  fi;
  for iVert in [1..Length(ListVert1)]
  do
    Add(eMat1, ListVert1[iVert]);
    jVert:=OnPoints(iVert, eGen);
    Add(eMat2, ListVert2[jVert]);
  od;
  test1:=Set(List(ListVert1, x->x[1]))=[0];
  test2:=Set(List(ListVert2, x->x[1]))=[0];
  if test1<>test2 then
    Error("We have a mixup here in the nature of the cones");
  fi;
  if test1=true then
    Add(eMat1, FacetOfInfinity(Length(ListVert1[1])));
    Add(eMat2, FacetOfInfinity(Length(ListVert1[1])));
  fi;
  rnk:=RankMat(eMat2);
  if rnk<>Length(ListVert1[1]) then
    Error("Rank Error");
  fi;
  Irr1:=[];
  Irr2:=[];
  rank:=0;
  pos:=1;
  while (rank<rnk)
  do
    TE:=ShallowCopy(Irr1);
    Add(TE, eMat1[pos]);
    if (RankMat(TE) > rank) then
      Irr1:=ShallowCopy(TE);
      Add(Irr2, eMat2[pos]);
      rank:=rank+1;
    fi;
    pos:=pos+1;
  od;
  TheMat:=Inverse(Irr1)*Irr2;
  if eMat1*TheMat<>eMat2 then
    return fail;
  fi;
  return TheMat;
end;


PersoGroup:=function(ListGen, eId)
  if Length(ListGen)=0 then
      return Group([eId]);
  else
      return Group(ListGen);
  fi;
end;


PersoGroupMatrix:=function(ListGen, n)
  if Length(ListGen)=0 then
    return Group([IdentityMat(n)]);
  else
    return Group(ListGen);
  fi;
end;


PersoGroupPerm:=function(ListGen)
  if Length(ListGen)=0 then
    return Group(());
  else
    return Group(ListGen);
  fi;
end;


SecondReduceGroupAction:=function(TheGroup, SmallSet)
  local ListGens, eGen, ListPos;
  ListGens:=[];
  if Length(SmallSet)<>Length(Set(SmallSet)) then
    Error("The input set has repetitions");
  fi;
  for eGen in GeneratorsOfGroup(TheGroup)
  do
    ListPos:=List(SmallSet, x->Position(SmallSet, OnPoints(x, eGen)));
    Add(ListGens, PermList(ListPos));
  od;
  return PersoGroupPerm(ListGens);
end;


# we complete a suitable vector family to
# a Z-spanning family of vectors
SubspaceCompletionInt:=function(eBasis, n)
  local RSE, d, A1, A2, TheProd, A1bis, i, j, FullBasis, TheCompletion;
  RSE:=SmithNormalFormIntegerMatTransforms(eBasis);
  d:=Length(eBasis);
  if Length(eBasis)=0 then
    return IdentityMat(n);
  fi;
  if RankMat(eBasis)<>d then
    Error("The vector family is not linearly independent, no way to complete");
  fi;
  A1:=RSE.rowtrans;
  A2:=RSE.coltrans;
  TheProd:=A1*eBasis*A2;
  for i in [1..d]
  do
    if TheProd[i][i]<>1 then
      Print("The basis B is independent but does not Z-span QB inter Z^n\n");
      Error("So we cannot extend it to a full Z-basis");
    fi;
  od;
  A1bis:=NullMat(n,n);
  for i in [1..d]
  do
    for j in [1..d]
    do
      A1bis[i][j]:=A1[i][j];
    od;
  od;
  for i in [d+1..n]
  do
    A1bis[i][i]:=1;
  od;
  FullBasis:=Inverse(A1bis)*Inverse(A2);
  TheCompletion:=FullBasis{[d+1..n]};
  return TheCompletion;
end;


# Find an orthonormal basis to eBasis
SubspaceCompletionRational:=function(eBasis, n)
  if Length(eBasis)=0 then
    return IdentityMat(n);
  fi;
  return NullspaceMat(TransposedMat(eBasis));
end;


__ProjectionFrame:=function(EXT)
  local rnk, EXTred;
  EXTred:=ShallowCopy(EXT);
  if TestConicness(EXTred) then
    Add(EXTred, FacetOfInfinity(Length(EXTred[1])) );
  fi;
  return ColumnReduction(EXTred, RankMat(EXTred)).Select;
end;


GetListPermGens:=function(ListVect, ListMatrGens)
  local eSperm, ListPermGens, ListVectSet, ListVectRed, eMatrGen, ListVectImg, eNewPerm;
  eSperm:=SortingPerm(ListVect);
  ListPermGens:=[];
  ListVectRed:=ColumnReduction(ListVect).EXT;
  ListVectSet:=Set(ListVect);
  for eMatrGen in ListMatrGens
  do
    ListVectImg:=ListVectSet*eMatrGen;
    if Set(ListVectImg)<>ListVectSet then
      Error("The matrix eMatrGen does not preserve ListVect");
    fi;
    eNewPerm:=eSperm*SortingPerm(ListVectImg)*Inverse(eSperm);
    Add(ListPermGens, eNewPerm);
    if FindTransformation(ListVectRed, ListVectRed, eNewPerm)=fail then
      Error("Bug in the code of GetListPermGens");
    fi;
  od;
  return ListPermGens;
end;


SYMPOL_PrintGroupStream:=function(output, n, GRP)
  local ListGen, eGen, i, j;
  ListGen:=GeneratorsOfGroup(GRP);
  AppendTo(output, n, " ", Length(ListGen), "\n");
  for eGen in ListGen
  do
    for i in [1..n]
    do
      j:=OnPoints(i, eGen);
      if j>n then
          Error("We have j=", j, " but n=", n);
      fi;
      AppendTo(output, " ", j-1);
    od;
    AppendTo(output, "\n");
  od;
end;


SYMPOL_PrintGroup:=function(eFile, n, GRP)
    local output;
    RemoveFileIfExist(eFile);
    output:=OutputTextFile(eFile, true);
    SYMPOL_PrintGroupStream(output, n, GRP);
    CloseStream(output);
end;


PersoRankMat:=function(ListVect)
  if Length(ListVect)=0 then
    return 0;
  fi;
  return RankMat(ListVect);
end;


GetTranslationClasses:=function(eMat)
  local n, ListGen, ListTransClass, ListStatus, FuncInsert, nbClass, IsFinished, iClass, eClass, eGen;
  n:=Length(eMat[1]);
  ListGen:=IdentityMat(n);
  ListTransClass:=[];
  ListStatus:=[];
  FuncInsert:=function(eVect)
    local fVect, eSol;
    for fVect in ListTransClass
    do
      eSol:=SolutionIntMat(eMat, eVect - fVect);
      if eSol<>fail then
        return;
      fi;
    od;
    Add(ListTransClass, eVect);
    Add(ListStatus, 0);
  end;
  FuncInsert(ListWithIdenticalEntries(n,0));
  while(true)
  do
    nbClass:=Length(ListTransClass);
    IsFinished:=true;
    for iClass in [1..nbClass]
    do
      if ListStatus[iClass]=0 then
        ListStatus[iClass]:=1;
	IsFinished:=false;
        eClass:=ListTransClass[iClass];
	for eGen in ListGen
	do
	  FuncInsert(eClass + eGen);
	od;
      fi;
    od;
    if IsFinished=true then
      break;
    fi;
    Print("Now |ListTransClass|=", Length(ListTransClass), "\n");
  od;
  return ListTransClass;
end;


# in this new version all incidences are encoded by sets.
__ProjectionLiftingFramework_Rational:=function(EXT, OneInc)
  local eSub, EXTproj, EXTproj2, FuncLift, RecReturn;
  eSub:=__ProjectionFrame(EXT);
  EXTproj:=List(EXT, x->x{eSub});
  EXTproj2:=EXTproj{OneInc};
  # TheInc is a subset of [1..Length(OneInc)]
  RecReturn:=rec(EXT:=EXT, OneInc:=OneInc);
  if RankMat(EXT)=2 then
    FuncLift:=function(TheInc)
      return Difference([1..2], OneInc);
    end;
    RecReturn.FuncLift:=FuncLift;
    return RecReturn;
  fi;
  FuncLift:=function(TheInc)
    local LINCR, VMA, ListCase, eExt, VO, EXT1, iD, EXT2, det12, iElt, EXTN, det1N, det2N, LV, test, S, iV, eV, eInc, VCE;
    LINCR:=EXTproj2{TheInc};
    if TestConicness(EXTproj) then
      Add(LINCR, FacetOfInfinity(Length(LINCR[1])) );
    fi;
    VMA:=NullspaceMat(TransposedMat(LINCR));
    if Length(VMA)<>2 then
      Error("We have an incoherence in ProjectionLiftingFramework");
    fi;
    ListCase:=[];
    for eExt in EXTproj
    do
      VO:=[VMA[1]*eExt, VMA[2]*eExt];
      if VO<>[0,0] then
        Add(ListCase, VO);
      fi;
    od;
    EXT1:=ListCase[1];
    iD:=2;
    while(true)
    do
      EXT2:=ListCase[iD];
      det12:=EXT2[2]*EXT1[1]-EXT2[1]*EXT1[2];
      if (det12<>0) then
        break;
      fi;
      iD:=iD+1;
    od;
    for iElt in [iD+1..Length(ListCase)]
    do
      EXTN:=ListCase[iElt];
      det1N:=EXTN[2]*EXT1[1]-EXTN[1]*EXT1[2];
      det2N:=EXTN[2]*EXT2[1]-EXTN[1]*EXT2[2];
      if (det1N*det2N>0) then
        if (det12*det1N>0) then
          EXT2:=ShallowCopy(EXTN);
          det12:=det1N;
        else
          EXT1:=ShallowCopy(EXTN);
          det12:=-det2N;
        fi;
      fi;
    od;
    if (det12>0) then
      LV:=[[-EXT1[2], EXT1[1]], [EXT2[2], -EXT2[1]]];
    else
      LV:=[[EXT1[2], -EXT1[1]], [-EXT2[2], EXT2[1]]];
    fi;
    test:=[1,1];
    S:=[];
    for iV in [1,2]
    do
      eV:=LV[iV];
      S[iV]:=eV[1]*VMA[1]+eV[2]*VMA[2];
      for eInc in EXTproj2
      do
        if S[iV]*eInc<>0 then
          test[iV]:=0;
        fi;
      od;
    od;
    VCE:=S[Position(test, 0)];
    return Filtered([1..Length(EXTproj)], x->EXTproj[x]*VCE=0);
  end;
  RecReturn.FuncLift:=FuncLift;
  return RecReturn;
end;


__ProjectionLiftingFramework:=function(EXT, OneInc)
  local Nval;
  if IsMatrixRational(EXT)=true then
    return __ProjectionLiftingFramework_Rational(EXT, OneInc);
  fi;
  Error("You have to build your own arithmetic");
end;


IntSqrt:=function(eVal)
    local LFact, LColl, TheRet, eColl, eMult, pVal, i;
    if eVal < 0 then
        return fail;
    fi;
    if eVal=1 then
        return 1;
    fi;
    LFact:=FactorsInt(eVal);
    LColl:=Collected(LFact);
    TheRet:=1;
    for eColl in LColl
    do
        eMult:=eColl[2] / 2;
        pVal:=eColl[1];
        if IsInt(eMult)=false then
            return fail;
        fi;
        for i in [1..eMult]
        do
            TheRet:=TheRet * pVal;
        od;
    od;
    return TheRet;
end;


RatSqrt:=function(eVal)
    local valDen, valNum;
    if eVal < 0 then
        return fail;
    fi;
    valDen := IntSqrt(DenominatorRat(eVal));
    valNum := IntSqrt(NumeratorRat(eVal));
    if valNum=fail or valDen=fail then
        return fail;
    fi;
    return valNum / valDen;
end;


LLLbasisReduction:=function(Latt)
    local GramMat, eRec, Pmat, LattRed;
    GramMat:=Latt * TransposedMat(Latt);
    eRec:=LLLReducedGramMat(GramMat);
    Pmat:=eRec.transformation;
    LattRed:=Pmat * Latt;
    return rec(LattRed:=LattRed, Pmat:=Pmat);
end;


# We want to consider the equation X A = B
# The equation is potentially underdefined and B in Z^*
# We look for the number d>0 such that for all B in Z^*
# the equation has a solution in Z^* / d
GetDenominatorQuotientSolution:=function(Amat)
    local AmatRed1, AmatRed2;
    AmatRed1:=GetZbasis(Amat);
    AmatRed2:=GetZbasisColumn(AmatRed1);
    return AbsInt(DeterminantMat(AmatRed2));
end;


#
# Add multiples of eVect so as 
EliminateSuperfluousPrimeDenominators:=function(eVect, ListVect)
    local dim, ListVect1, ListVect2, ListVect3, ListVect4, TheCompl, TheBasis, eSol, eSolRed, TheRet;
    dim:=Length(eVect);
    if Length(ListVect)=0 then
        # nothing possible really.
        return eVect;
    fi;
    if RankMat(ListVect) = Length(eVect) then
        return ListWithIdenticalEntries(dim,0);
    fi;
    ListVect1:=List(ListVect, RemoveFraction);
    ListVect2:=NullspaceMat(TransposedMat(ListVect1));
    ListVect3:=List(ListVect2, RemoveFraction);
    ListVect4:=NullspaceIntMat(TransposedMat(ListVect3));
    #
    TheCompl:=SubspaceCompletionInt(ListVect4, dim);
    TheBasis:=Concatenation(TheCompl, ListVect4);
    eSol:=eVect * Inverse(TheBasis);
    eSolRed:=eSol{[1..Length(TheCompl)]};
    TheRet:=eSolRed * TheCompl;
    return TheRet;
end;


EliminateSuperfluousPrimeDenominators_Matrix:=function(eMat, ListMat)
    local eVect, ListVect, TheVect, TheMat;
    eVect:=MatrixToVector(eMat);
    ListVect:=List(ListMat, MatrixToVector);
    TheVect:=EliminateSuperfluousPrimeDenominators(eVect, ListVect);
    TheMat:=VectorToMatrix(TheVect);
    return TheMat;
end;
