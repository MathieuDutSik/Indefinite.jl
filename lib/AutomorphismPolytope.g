FileGRP_ComputeAut_ListMat_Subset_EXT:=Filename(DirectoriesPackagePrograms("indefinite"),"GRP_ListMat_Subset_EXT_Automorphism");
FileGRP_TestEquivalence_ListMat_Subset_EXT:=Filename(DirectoriesPackagePrograms("indefinite"),"GRP_ListMat_Subset_EXT_Isomorphism");
FileGRP_Invariant_ListMat_Subset_EXT:=Filename(DirectoriesPackagePrograms("indefinite"),"GRP_ListMat_Subset_EXT_Invariant");


ShortVectorDutourVersion:=function(GramMat, eNorm)
  local H, eNormRed, GramMatRed, SHV, ListVector, iV, eV;
  H:=RemoveFractionMatrixPlusCoef(GramMat);
  eNormRed:=LowerInteger(eNorm*H.TheMult);
  GramMatRed:=H.TheMat;
  if IsPositiveDefiniteSymmetricMatrix(GramMatRed)=false then
    Error("Matrix should be positive definite");
  fi;
  SHV:=ShortestVectors(GramMatRed, eNormRed);
  ListVector:=[];
  for iV in [1..Length(SHV.vectors)]
  do
    if SHV.norms[iV] <= eNormRed then
      eV:=SHV.vectors[iV];
      Add(ListVector, eV);
      Add(ListVector, -eV);
    fi;
  od;
  return ListVector;
end;


__VectorConfigurationFullDim_ScalarMat_AddMat:=function(EXT, ListAddMat)
  local n, Qmat, eEXT, Qinv, ScalarMat, fEXT, eLine, ListMat, LVal;
  n:=Length(EXT[1]);
  if Length(EXT)<>Length(Set(EXT)) then
    Error("The vertex list has some repetition");
  fi;
  if RankMat(EXT)<>n then
    Error("The polytope is not of ful rank");
  fi;
  Qmat:=NullMat(n,n);
  for eEXT in EXT
  do
    Qmat:=Qmat+TransposedMat([eEXT])*[eEXT];
  od;
  Qinv:=Inverse(Qmat);
  ListMat:=Concatenation([Qinv], ListAddMat);
  ScalarMat:=[];
  for eEXT in EXT
  do
    eLine:=[];
    for fEXT in EXT
    do
      LVal:=List(ListMat, x->eEXT*x*fEXT);
      Add(eLine, LVal);
    od;
    Add(ScalarMat, eLine);
  od;
  return ScalarMat;
end;


__VectorConfiguration_Invariant_GetTools:=function(EXT, TheLimit)
  local eRec, EXTred, eSelect, n, Qmat, eEXT, Qinv, ListDiagVal, PreListOffDiag, ListOffDiag, iVert, nbVert, eProd;
  eRec:=ColumnReduction(EXT, RankMat(EXT));
  EXTred:=eRec.EXT;
  eSelect:=eRec.Select;
  n:=Length(EXTred[1]);
  Qmat:=NullMat(n,n);
  for eEXT in EXTred
  do
    Qmat:=Qmat+TransposedMat([eEXT])*[eEXT];
  od;
  Qinv:=Inverse(Qmat);
  nbVert:=Length(EXT);
  ListDiagVal:=Set(List([1..nbVert], x->EXTred[x]*Qinv*EXTred[x]));
  if nbVert <= TheLimit then
    PreListOffDiag:=[];
    for iVert in [1..nbVert-1]
    do
      eProd:=EXTred[iVert]*Qinv;
      Append(PreListOffDiag, List([iVert+1..nbVert], x->eProd*EXTred[x]));
    od;
    ListOffDiag:=Set(PreListOffDiag);
  else
    ListOffDiag:=[];
  fi;
  return rec(eSelect:=eSelect, Qinv:=Qinv, TheLimit:=TheLimit, EXTred:=EXTred,
             ListDiagVal:=ListDiagVal,
             ListOffDiag:=ListOffDiag);
end;


__VectorConfiguration_Invariant_Compute:=function(eTool, EXT)
  local EXTred, n, Qmat, eEXT, Qinv, ListValDiag, ListNbDiag, ListValOff, ListNbOff, iVert, jVert, eScal, pos, nbVert, PreListValDiag, PreListNbDiag, PreListValOff, PreListNbOff, ePerm;
  EXTred:=List(EXT, x->x{eTool.eSelect});
  n:=Length(EXTred[1]);
  nbVert:=Length(EXT);
  PreListValDiag:=[];
  PreListNbDiag:=[];
  for iVert in [1..nbVert]
  do
    eScal:=EXTred[iVert]*eTool.Qinv*EXTred[iVert];
    pos:=Position(PreListValDiag, eScal);
    if pos=fail then
      Add(PreListValDiag, eScal);
      Add(PreListNbDiag, 1);
    else
      PreListNbDiag[pos]:=PreListNbDiag[pos]+1;
    fi;
  od;
  ePerm:=SortingPerm(PreListValDiag);
  ListValDiag:=Permuted(PreListValDiag, ePerm);
  ListNbDiag:=Permuted(PreListNbDiag, ePerm);
  PreListValOff:=[];
  PreListNbOff:=[];
  if nbVert< eTool.TheLimit then
    for iVert in [1..nbVert-1]
    do
      for jVert in [iVert+1..nbVert]
      do
        eScal:=EXTred[iVert]*eTool.Qinv*EXTred[jVert];
        pos:=Position(PreListValOff, eScal);
        if pos=fail then
          Add(PreListValOff, eScal);
          Add(PreListNbOff, 1);
        else
          PreListNbOff[pos]:=PreListNbOff[pos]+1;
        fi;
      od;
    od;
  fi;
  ePerm:=SortingPerm(PreListValOff);
  ListValOff:=Permuted(PreListValOff, ePerm);
  ListNbOff:=Permuted(PreListNbOff, ePerm);
  return rec(n:=n,
             ListValDiag:=ListValDiag, ListNbDiag:=ListNbDiag,
             ListValOff:=ListValOff, ListNbOff:=ListNbOff);
end;


__VectorConfiguration_Invariant_ComputeAdvanced:=function(eTool, eInc)
  local nbVert, ListScal, nbDiag, nbOff, eVect1, eVect2, eVect3, diffInc, eVal, eV, eScal, pos, eLen, diffLen, eProd, i, j;
  nbVert:=Length(eTool.EXTred);
  nbDiag:=Length(eTool.ListDiagVal);
  eVect1:=ListWithIdenticalEntries(nbDiag,0);
  for eVal in eInc
  do
    eV:=eTool.EXTred[eVal];
    eScal:=eV * eTool.Qinv * eV;
    pos:=Position(eTool.ListDiagVal, eScal);
    eVect1[pos]:=eVect1[pos]+1;
  od;
  if eTool.TheLimit < nbVert then
    return eVect1;
  fi;
  nbOff :=Length(eTool.ListOffDiag);
  eVect2:=ListWithIdenticalEntries(nbOff ,0);
  eVect3:=ListWithIdenticalEntries(nbOff ,0);
  diffInc:=Difference([1..nbVert], eInc);
  eLen:=Length(eInc);
  diffLen:=Length(diffInc);
  for i in [1..eLen]
  do
    eProd:=eTool.EXTred[eInc[i]] * eTool.Qinv;
    for j in [i+1..eLen]
    do
      eScal:=eProd * eTool.EXTred[eInc[j]];
      pos:=Position(eTool.ListOffDiag, eScal);
      eVect2[pos]:=eVect2[pos]+1;
    od;
  od;
  for i in [1..eLen]
  do
    eProd:=eTool.EXTred[eInc[i]] * eTool.Qinv;
    for j in [1..diffLen]
    do
      eScal:=eProd * eTool.EXTred[diffInc[j]];
      pos:=Position(eTool.ListOffDiag, eScal);
      eVect3[pos]:=eVect3[pos]+1;
    od;
  od;
  return Concatenation(eVect1, eVect2, eVect3);
end;


__VectorConfiguration_Invariant:=function(EXT, TheLimit)
  local eTool;
  eTool:=__VectorConfiguration_Invariant_GetTools(EXT, TheLimit);
  return __VectorConfiguration_Invariant_Compute(eTool, EXT);
end;


LinPolytope_Invariant:=function(EXT)
  local TheLimit;
  TheLimit:=500;
  return __VectorConfiguration_Invariant(EXT, TheLimit);
end;


Get_RecScalColor:=function(EXT, GramMat)
  local GetScalarColor, GetLineColor, nbVert;
  GetScalarColor:=function(i,j)
    return EXT[i]*GramMat*EXT[j];
  end;
  GetLineColor:=function(i)
    local eVect;
    eVect:=EXT[i]*GramMat;
    return EXT*eVect;
  end;
  nbVert:=Length(EXT);
  return rec(n:=nbVert,
             GetScalarColor:=GetScalarColor,
             GetLineColor:=GetLineColor);
end;


LinPolytope_Automorphism_Scalable:=function(EXT, GramMat)
  local eRecScalColor;
  eRecScalColor:=Get_RecScalColor(EXT, GramMat);
  return AutomorphismGroupColoredGraph_Scalable(eRecScalColor);
end;


LinPolytope_Automorphism_Simple:=function(EXT, GramMat)
  local ScalarMat;
  ScalarMat:=EXT*GramMat*TransposedMat(EXT);
  return AutomorphismGroupColoredGraph(ScalarMat);
end;


LinPolytope_Automorphism_GramMat:=function(EXT, GramMat)
  if Length(EXT)<700 then
    return LinPolytope_Automorphism_Simple(EXT, GramMat);
  fi;
  return LinPolytope_Automorphism_Scalable(EXT, GramMat);
end;


Get_QinvMatrix:=function(EXT)
  local n, Qmat, eEXT;
  n:=Length(EXT[1]);
  Qmat:=NullMat(n,n);
  for eEXT in EXT
  do
    Qmat:=Qmat+TransposedMat([eEXT])*[eEXT];
  od;
  return RemoveFractionMatrix(Inverse(Qmat));
end;


LinPolytope_Automorphism:=function(EXT)
  local EXTred, n, Qmat, eEXT, Qinv;
  EXTred:=ColumnReduction(EXT).EXT;
  Qinv:=Get_QinvMatrix(EXTred);
  return LinPolytope_Automorphism_GramMat(EXTred, Qinv);
end;


Get_RecScalColor_Subset_AddMat:=function(EXT, EXTsub, ListAddMat)
  local eSet, GetScalarColor, GetLineColor, nbVert;
  eSet:=Set(List(EXTsub, x->Position(EXT, x)));
  nbVert:=Length(EXT);
  GetScalarColor:=function(i, j)
    local eLine, eMat;
    eLine:=[];
    for eMat in ListAddMat
    do
      Add(eLine, EXT[i] * eMat * EXT[j]);
    od;
    if i=j then
      if i in eSet then
        Add(eLine, 1);
      else
        Add(eLine, 0);
      fi;
    fi;
    return eLine;
  end;
  GetLineColor:=function(i)
    return List([1..nbVert], x->GetScalarColor(x, i));
  end;
  return rec(n:=nbVert,
             GetScalarColor:=GetScalarColor,
             GetLineColor:=GetLineColor);
end;


WriteData_ListMat_Subset:=function(output, EXT, EXTsub, ListMat)
    local ListVal, eEXT, eVal, eMat, i;
    ListVal:=[];
    for eEXT in EXT
    do
        eVal:=0;
        if Position(EXTsub, eEXT)<>fail then
            eVal:=1;
        fi;
        Add(ListVal, eVal);
    od;
    AppendTo(output, Length(ListMat), "\n");
    for eMat in ListMat
    do
        CPP_WriteMatrix(output, eMat);
    od;
    CPP_WriteMatrix(output, EXT);
    AppendTo(output, Length(EXT), "\n");
    for i in [1..Length(EXT)]
    do
        AppendTo(output, " ", ListVal[i]);
    od;
    AppendTo(output, "\n");
end;


GetScalarMatrix_PolytopeStabSubset_AddMat:=function(EXT, EXTsub, ListAddMat)
  local eSet, EXTred, ScalarMat, nbVert, RedoneScalarMat, eLine, iVert, jVert, eValMatr, eVal;
  eSet:=Set(List(EXTsub, x->Position(EXT, x)));
  if RankMat(EXT)<>Length(EXT[1]) then
    Error("Rank error in _AddMat function");
  fi;
  ScalarMat:=__VectorConfigurationFullDim_ScalarMat_AddMat(EXT, ListAddMat);
  nbVert:=Length(EXT);
  RedoneScalarMat:=[];
  for iVert in [1..nbVert]
  do
    eLine:=[];
    for jVert in [1..nbVert]
    do
      eValMatr:=ScalarMat[iVert][jVert];
      if iVert<>jVert then
        Add(eLine, eValMatr);
      else
        if iVert in eSet then
          eVal:=1;
        else
          eVal:=0;
        fi;
        Add(eLine, [eValMatr, eVal]);
      fi;
    od;
    Add(RedoneScalarMat, eLine);
  od;
  return RedoneScalarMat;
end;


GetScalarMatrixInvariant_PolytopeStabSubset_AddMat_GAP:=function(EXT, EXTsub, ListAddMat)
  local eSet, nbVert, GetValue, PreListValues, ListValues, iVert, jVert, eVal, pos, nbValues, ListOccDiag, ListOccOff;
  eSet:=Set(List(EXTsub, x->Position(EXT, x)));
  if RankMat(EXT)<>Length(EXT[1]) then
    Error("Rank error in _AddMat function");
  fi;
  nbVert:=Length(EXT);
  GetValue:=function(i, j)
    local eLine, eMat;
    eLine:=[];
    for eMat in ListAddMat
    do
      Add(eLine, EXT[i] * eMat * EXT[j]);
    od;
    if i=j then
      if i in eSet then
        Add(eLine, 1);
      else
        Add(eLine, 0);
      fi;
    fi;
    return eLine;
  end;
  PreListValues:=[];
  for iVert in [1..nbVert]
  do
    for jVert in [1..nbVert]
    do
      eVal:=GetValue(iVert, jVert);
      if Position(PreListValues, eVal)=fail then
        Add(PreListValues, eVal);
      fi;
    od;
  od;
  ListValues:=Set(PreListValues);
  nbValues:=Length(ListValues);
  ListOccDiag:=ListWithIdenticalEntries(nbValues, 0);
  for iVert in [1..nbVert]
  do
    eVal:=GetValue(iVert, iVert);
    pos:=Position(ListValues, eVal);
    ListOccDiag[pos]:=ListOccDiag[pos] + 1;
  od;
  ListOccOff:=ListWithIdenticalEntries(nbValues, 0);
  for iVert in [1..nbVert]
  do
    for jVert in [1..nbVert]
    do
      if iVert<>jVert then
        eVal:=GetValue(iVert, jVert);
        pos:=Position(ListValues, eVal);
        ListOccOff[pos]:=ListOccOff[pos] + 1;
      fi;
    od;
  od;
  return rec(ListValues:=ListValues, ListOccDiag:=ListOccDiag, ListOccOff:=ListOccOff);
end;


GetScalarMatrixInvariant_PolytopeStabSubset_AddMat:=function(EXT, EXTsub, ListAddMat)
    local ListMat, FileI, FileO, output, eInv;
    ListMat:=Concatenation([Get_QinvMatrix(EXT)], ListAddMat);
    FileI:=Filename(POLYHEDRAL_tmpdir,"ListMatSubset_auto.in");
    FileO:=Filename(POLYHEDRAL_tmpdir,"ListMatSubset_auto.out");
    output:=OutputTextFile(FileI, true);
    WriteData_ListMat_Subset(output, EXT, EXTsub, ListMat);
    CloseStream(output);
    #
    Exec(FileGRP_Invariant_ListMat_Subset_EXT, " ", FileI, " ", FileO);
    eInv:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    return eInv;
end;


LinPolytope_AutomorphismStabSubset_AddMat:=function(EXT, EXTsub, ListAddMat)
    local output, ListMat, GRP, FileI, FileO;
    ListMat:=Concatenation([Get_QinvMatrix(EXT)], ListAddMat);
    FileI:=Filename(POLYHEDRAL_tmpdir,"ListMatSubset_auto.in");
    FileO:=Filename(POLYHEDRAL_tmpdir,"ListMatSubset_auto.out");
    output:=OutputTextFile(FileI, true);
    WriteData_ListMat_Subset(output, EXT, EXTsub, ListMat);
    CloseStream(output);
    #
    Exec(FileGRP_ComputeAut_ListMat_Subset_EXT, " ", FileI, " ", FileO);
    GRP:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    return GRP;
end;


LinPolytope_IsomorphismStabSubset_AddMat:=function(EXT1, EXTsub1, EXT2, EXTsub2, ListAddMat1, ListAddMat2)
    local RedoneScalarMat1, RedoneScalarMat2, eEquiv, FileI, FileO, ListMat1, ListMat2, output, TheOut, UseMethod;
    FileI:=Filename(POLYHEDRAL_tmpdir,"ListMatSubset_isom.in");
    FileO:=Filename(POLYHEDRAL_tmpdir,"ListMatSubset_isom.out");
    ListMat1:=Concatenation([Get_QinvMatrix(EXT1)], ListAddMat1);
    ListMat2:=Concatenation([Get_QinvMatrix(EXT2)], ListAddMat2);
    output:=OutputTextFile(FileI, true);
    WriteData_ListMat_Subset(output, EXT1, EXTsub1, ListMat1);
    WriteData_ListMat_Subset(output, EXT2, EXTsub2, ListMat2);
    CloseStream(output);
    Exec(FileGRP_TestEquivalence_ListMat_Subset_EXT, " ", FileI, " ", FileO);
    TheOut:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    if TheOut=false then
        return false;
    else
        return PermList(TheOut);
    fi;
end;


LinPolytope_Isomorphism_Simple:=function(EXT1, GramMat1, EXT2, GramMat2)
  local eEquiv, ScalarMat1, ScalarMat2;
  ScalarMat1:=EXT1*GramMat1*TransposedMat(EXT1);
  ScalarMat2:=EXT2*GramMat2*TransposedMat(EXT2);
  eEquiv:=IsIsomorphicColoredGraph(ScalarMat1, ScalarMat2);
  if eEquiv=false then
    return false;
  fi;
  return PermList(eEquiv);
end;


LinPolytope_Isomorphism_Scalable:=function(EXT1, GramMat1, EXT2, GramMat2)
  local eRecScalColor1, eRecScalColor2, eEquiv;
  eRecScalColor1:=Get_RecScalColor(EXT1, GramMat1);
  eRecScalColor2:=Get_RecScalColor(EXT2, GramMat2);
  eEquiv:=IsIsomorphicColoredGraph_Scalable(eRecScalColor1, eRecScalColor2);
  if eEquiv=false then
    return false;
  fi;
  return PermList(eEquiv);
end;


LinPolytope_Isomorphism_GramMat:=function(EXT1, GramMat1, EXT2, GramMat2)
  if Length(EXT1)<>Length(EXT2) then
    return false;
  fi;
  if Length(EXT1)<700 then
    return LinPolytope_Isomorphism_Simple(EXT1, GramMat1, EXT2, GramMat2);
  fi;
  return LinPolytope_Isomorphism_Scalable(EXT1, GramMat1, EXT2, GramMat2);
end;


LinPolytope_Isomorphism:=function(EXT1, EXT2)
  local EXTred1, EXTred2, Qinv1, Qinv2;
  EXTred1:=ColumnReduction(EXT1).EXT;
  EXTred2:=ColumnReduction(EXT2).EXT;
  Qinv1:=Get_QinvMatrix(EXTred1);
  Qinv2:=Get_QinvMatrix(EXTred2);
  return LinPolytope_Isomorphism_GramMat(EXTred1, Qinv1, EXTred2, Qinv2);
end;


LinPolytopeIntegral_Isomorphism_Exhaustive:=function(EXT1, EXT2)
  local eEquiv, GRP, eElt, nEquiv, eBigMat;
  eEquiv:=LinPolytope_Isomorphism(EXT1, EXT2);
  if eEquiv=false then
    return false;
  fi;
  GRP:=LinPolytope_Automorphism(EXT1);
  for eElt in GRP
  do
    nEquiv:=eElt*eEquiv;
    eBigMat:=FindTransformation(EXT1, EXT2, nEquiv);
    if IsIntegralMat(eBigMat) then
      return eBigMat;
    fi;
  od;
  return false;
end;


#
# Function below is buggy. It does not compute correctly the
# Delaunay in dimension 5 from the Delaunay in dimension 6.
KernelLinPolytopeIntegral_Isomorphism_Subspaces:=function(EXT1, EXT2, GRP2, eEquiv)
  local n, eBasis1, eBasis2, EXTbas1, EXTbas2, TheMatEquiv, ListMatrGens, eGen, TheMat, GRPspace, eLatt1, eLatt2, eRec1, eRec2, eSpaceEquiv, eMatFinal;
  Print("Begin KernelLinPolytopeIntegral_Isomorphism_Subspaces\n");
  n:=Length(EXT1[1]);
  eBasis1:=GetZbasis(EXT1);
  eBasis2:=GetZbasis(EXT2);
  EXTbas1:=List(EXT1, x->SolutionMat(eBasis1, x));
  EXTbas2:=List(EXT2, x->SolutionMat(eBasis2, x));
  TheMatEquiv:=FindTransformation(EXTbas1, EXTbas2, eEquiv);
  Print("After FindTransformation\n");
  ListMatrGens:=[];
  for eGen in GeneratorsOfGroup(GRP2)
  do
    TheMat:=FindTransformation(EXTbas2, EXTbas2, eGen);
    Add(ListMatrGens, TheMat);
  od;
  if Length(ListMatrGens)=0 then
      GRPspace:=Group([IdentityMat(n)]);
  else
      GRPspace:=Group(ListMatrGens);
  fi;
  eLatt1:=Inverse(eBasis1)*TheMatEquiv;
  eLatt2:=Inverse(eBasis2);
  eRec1:=RemoveFractionMatrixPlusCoef(eLatt1);
  eRec2:=RemoveFractionMatrixPlusCoef(eLatt2);
  if eRec1.TheMult<>eRec2.TheMult then
    return false;
  fi;
#  Print("Before call to LinearSpace_Equivalence\n");
  eSpaceEquiv:=LinearSpace_Equivalence(GRPspace, eRec1.TheMat, eRec2.TheMat);
#  Print("After call to LinearSpace_Equivalence\n");
  if eSpaceEquiv=fail then
    return false;
  fi;
  eMatFinal:=Inverse(eBasis1)*TheMatEquiv*eSpaceEquiv*eBasis2;
  if IsIntegralMat(eMatFinal)=false then
    Error("eMatFinal is not integral, BUG");
  fi;
  if Set(EXT1*eMatFinal)<>Set(EXT2) then
    Error("eMatFinal does not map the polytopes, BUG");
  fi;
  return eMatFinal;
end;


LinPolytopeIntegral_Isomorphism_Subspaces:=function(EXT1, EXT2)
  local n, eEquiv, GRP2;
  n:=Length(EXT1[1]);
  if Length(EXT2[1])<>n then
    Error("The dimension of EXT1 and EXT2 are not the same");
  fi;
  if RankMat(EXT1)<>n then
    Error("EXT1 is not full dimensional");
  fi;
  if RankMat(EXT2)<>n then
    Error("EXT2 is not full dimensional");
  fi;
  eEquiv:=LinPolytope_Isomorphism(EXT1, EXT2);
  if eEquiv=false then
    return false;
  fi;
  GRP2:=LinPolytope_Automorphism(EXT2);
  Print("|GRP2|=", Order(GRP2), "\n");
  return KernelLinPolytopeIntegral_Isomorphism_Subspaces(EXT1, EXT2, GRP2, eEquiv);
end;


KernelLinPolytopeIntegral_Automorphism_Subspaces:=function(EXT, GRP)
  local n, eBasis, EXTbas, ListPermGens, ListMatrGens, eGen, TheMat, GRPmatr, LattToStab, eStab, eList, ePerm, eMatr, GRPpermFinal, GRPmatrFinal;
  n:=Length(EXT[1]);
  if RankMat(EXT)<>n then
    Error("Wrong rank for LinPolytopeIntegral_Automorphism");
  fi;
  eBasis:=GetZbasis(EXT);
  EXTbas:=List(EXT, x->SolutionMat(eBasis, x));
#  Print("|GRP|=", Order(GRP), "\n");
  ListPermGens:=GeneratorsOfGroup(GRP);
  ListMatrGens:=[];
  for eGen in ListPermGens
  do
    TheMat:=FindTransformation(EXTbas, EXTbas, eGen);
    Add(ListMatrGens, TheMat);
  od;
  if Length(ListMatrGens)=0 then
      GRPmatr:=Group([IdentityMat(n)]);
  else
      GRPmatr:=Group(ListMatrGens);
  fi;
  LattToStab:=RemoveFractionMatrix(Inverse(eBasis));
  #
  eStab:=LinearSpace_Stabilizer(GRPmatr, LattToStab);
  ListPermGens:=[];
  ListMatrGens:=[];
  for eGen in GeneratorsOfGroup(eStab)
  do
    eList:=List(EXTbas, x->Position(EXTbas, x*eGen));
    ePerm:=PermList(eList);
    eMatr:=FindTransformation(EXT, EXT, ePerm);
    if IsIntegralMat(eMatr)=false then
      Error("Some bugs to resolve");
    fi;
    Add(ListPermGens, ePerm);
    Add(ListMatrGens, eMatr);
  od;
  if Length(ListMatrGens)=0 then
      GRPmatrFinal:=Group([IdentityMat(n)]);
  else
      GRPmatrFinal:=Group(ListMatrGens);
  fi;
  GRPpermFinal:=PersoGroupPerm(ListPermGens);
  SetSize(GRPmatrFinal, Order(GRPpermFinal));
  return rec(GRPperm:=GRPpermFinal, GRPmatr:=GRPmatrFinal);
end;


LinPolytopeIntegral_Isomorphism:=function(EXT1, EXT2)
  local GRP1;
  GRP1:=LinPolytope_Automorphism(EXT1);
  if Order(GRP1)<=10000 then
    return LinPolytopeIntegral_Isomorphism_Exhaustive(EXT1, EXT2);
  fi;
  return LinPolytopeIntegral_Isomorphism_Subspaces(EXT1, EXT2);
end;
