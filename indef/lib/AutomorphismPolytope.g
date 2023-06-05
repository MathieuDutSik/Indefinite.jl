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


LinPolytope_Automorphism_GramMat:=function(EXT, GramMat)
    local EXT_oscar, GramMat_oscar;
    EXT_oscar:=MatrixToOscar(EXT);
    GramMat_oscar:=MatrixToOscar(GramMat);
    return Julia.Indefinite.GRP_LinPolytope_Automorphism_GramMat(EXT_oscar, GramMat_oscar);
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


ParseMyOscarPermIsomorphism:=function(eEquiv_oscar)
    local TheMat;
    TheMat:=ReadOscarMatrix(eEquiv_oscar);
    if Length(TheMat) =0 then
        return false;
    else
        return PermList(TheMat[1]);
    fi;
end;



LinPolytope_Isomorphism_GramMat:=function(EXT1, GramMat1, EXT2, GramMat2)
    local EXT1_oscar, GramMat1_oscar, EXT2_oscar, GramMat2_oscar, eEquiv_oscar;
    EXT1_oscar:=MatrixToOscar(EXT1);
    GramMat1_oscar:=MatrixToOscar(GramMat1);
    EXT2_oscar:=MatrixToOscar(EXT2);
    GramMat2_oscar:=MatrixToOscar(GramMat2);
    eEquiv_oscar:=Julia.Indefinite.GRP_LinPolytope_Isomorphism_GramMat(EXT1_oscar, GramMat1_oscar, EXT2_oscar, GramMat2_oscar);
    return ParseMyOscarPermIsomorphism(eEquiv_oscar);
end;


LinPolytope_Isomorphism:=function(EXT1, EXT2)
  local EXTred1, EXTred2, Qinv1, Qinv2;
  EXTred1:=ColumnReduction(EXT1).EXT;
  EXTred2:=ColumnReduction(EXT2).EXT;
  Qinv1:=Get_QinvMatrix(EXTred1);
  Qinv2:=Get_QinvMatrix(EXTred2);
  return LinPolytope_Isomorphism_GramMat(EXTred1, Qinv1, EXTred2, Qinv2);
end;


GetScalarMatrixInvariant_Polytope_AddMat:=function(EXT, ListAddMat)
    local EXT_oscar, Qinv, ListMat, ListMat_oscar, Vdiag, Vdiag_oscar;
    EXT_oscar:=MatrixToOscar(EXT);
    Qinv:=Get_QinvMatrix(EXT);
    ListMat:=Concatenation([Get_QinvMatrix(EXT)], ListAddMat);
    ListMat_oscar:=ListMatrixToOscar(ListMat);
    Vdiag:=[ListWithIdenticalEntries(Length(EXT),1)];
    Vdiag_oscar:=MatrixToOscar(Vdiag);
    return Julia.Indefinite.GRP_ListMat_Subset_EXT_Invariant(EXT_oscar, ListMat_oscar, Vdiag_oscar);
end;


LinPolytope_Invariant:=function(EXT)
    local EXTred;
    EXTred:=ColumnReduction(EXT).EXT;
    return GetScalarMatrixInvariant_Polytope_AddMat(EXTred, []);
end;


LinPolytope_Automorphism_AddMat:=function(EXT, ListAddMat)
    local EXT_oscar, Qinv, ListMat, ListMat_oscar, Vdiag, Vdiag_oscar;
    EXT_oscar:=MatrixToOscar(EXT);
    Qinv:=Get_QinvMatrix(EXT);
    ListMat:=Concatenation([Get_QinvMatrix(EXT)], ListAddMat);
    ListMat_oscar:=ListMatrixToOscar(ListMat);
    Vdiag:=[ListWithIdenticalEntries(Length(EXT),1)];
    Vdiag_oscar:=MatrixToOscar(Vdiag);
    return Julia.Indefinite.GRP_ListMat_Subset_EXT_Automorphism(EXT_oscar, ListMat_oscar, Vdiag_oscar);
end;


LinPolytope_Isomorphism_AddMat:=function(EXT1, EXT2, ListAddMat1, ListAddMat2)
    local EXT1_oscar, Qinv1, ListMat1, ListMat1_oscar, EXT2_oscar, Qinv2, ListMat2, ListMat2_oscar, eEquiv_oscar, Vdiag1, Vdiag2, Vdiag1_oscar, Vdiag2_oscar;
    EXT1_oscar:=MatrixToOscar(EXT1);
    Qinv1:=Get_QinvMatrix(EXT1);
    ListMat1:=Concatenation([Get_QinvMatrix(EXT1)], ListAddMat1);
    ListMat1_oscar:=ListMatrixToOscar(ListMat1);
    EXT2_oscar:=MatrixToOscar(EXT2);
    Qinv2:=Get_QinvMatrix(EXT2);
    ListMat2:=Concatenation([Get_QinvMatrix(EXT2)], ListAddMat2);
    ListMat2_oscar:=ListMatrixToOscar(ListMat2);
    Vdiag1:=NullMat(1, Length(EXT1));
    Vdiag2:=NullMat(1, Length(EXT2));
    Vdiag1_oscar:=MatrixToOscar(Vdiag1);
    Vdiag2_oscar:=MatrixToOscar(Vdiag2);
    eEquiv_oscar:=Julia.Indefinite.GRP_ListMat_Subset_EXT_Isomorphism(EXT1_oscar, ListMat1_oscar, Vdiag1_oscar, EXT2_oscar, ListMat2_oscar, Vdiag2_oscar);
    return ParseMyOscarPermIsomorphism(eEquiv_oscar);
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
    if IndefinitePrint then
        Print("Begin KernelLinPolytopeIntegral_Isomorphism_Subspaces\n");
    fi;
    n:=Length(EXT1[1]);
    eBasis1:=GetZbasis(EXT1);
    eBasis2:=GetZbasis(EXT2);
    EXTbas1:=List(EXT1, x->SolutionMat(eBasis1, x));
    EXTbas2:=List(EXT2, x->SolutionMat(eBasis2, x));
    TheMatEquiv:=FindTransformation(EXTbas1, EXTbas2, eEquiv);
    if IndefinitePrint then
        Print("After FindTransformation\n");
    fi;
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
    if IndefinitePrint then
        Print("Before call to LinearSpace_Equivalence\n");
    fi;
    eSpaceEquiv:=LinearSpace_Equivalence(GRPspace, eRec1.TheMat, eRec2.TheMat);
    if IndefinitePrint then
        Print("After call to LinearSpace_Equivalence\n");
    fi;
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
    if IndefinitePrint then
        Print("|GRP2|=", Order(GRP2), "\n");
    fi;
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
  if IndefinitePrint then
      Print("|GRP|=", Order(GRP), "\n");
  fi;
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
