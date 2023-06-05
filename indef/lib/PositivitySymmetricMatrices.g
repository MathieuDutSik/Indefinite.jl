IsSymmetricMatrix:=function(TheMat)
    return TransposedMat(TheMat) = TheMat;
end;


SymmetricExtractSubMatrix:=function(SymMat, eSet)
  return List(SymMat{eSet}, x->x{eSet});
end;


IsPositiveDefiniteSymmetricMatrix:=function(SymMat)
  local iCol;
  if TransposedMat(SymMat)<>SymMat then
    Error("hum, the matrix should be symmetric?!");
  fi;
  for iCol in [1..Length(SymMat)]
  do
    if DeterminantMat(SymmetricExtractSubMatrix(SymMat, [1..iCol]))<=0 then
      return false;
    fi;
  od;
  return true;
end;


NullspaceReduction:=function(SymMat)
  local n, TransfMat, DimKern, RR, iCol, V, SymMat2, MatExt;
  if SymMat<>TransposedMat(SymMat) then
    Error("The matrix is not symmetric");
  fi;
  n:=Length(SymMat);
  TransfMat:=ShallowCopy(NullspaceMat(SymMat));
  DimKern:=Length(TransfMat);
  RR:=ColumnReduction(TransfMat, DimKern);
  for iCol in [1..n]
  do
    if Position(RR.Select, iCol)=fail then
      V:=ListWithIdenticalEntries(n, 0);
      V[iCol]:=1;
      Add(TransfMat, ShallowCopy(V));
    fi;
  od;
  SymMat2:=TransfMat*SymMat*TransposedMat(TransfMat);
  MatExt:=SymmetricExtractSubMatrix(SymMat2, [1+DimKern..n]);
  return rec(RedMat:=SymMat2, Transform:=TransfMat, NonDegenerate:=MatExt);
end;


IsPositiveSymmetricMatrix:=function(SymMat)
  return IsPositiveDefiniteSymmetricMatrix(NullspaceReduction(SymMat).NonDegenerate);
end;


# the matrix should be non-zero.
FindNonIsotropeVector:=function(SymMat)
  local n, i, j, VECT;
  n:=Length(SymMat);
  for i in [1..n]
  do
    if SymMat[i][i]<>0 then
      VECT:=ListWithIdenticalEntries(n, 0);
      VECT[i]:=1;
      return VECT;
    fi;
  od;
  for i in [1..n-1]
  do
    for j in [i+1..n]
    do
      if SymMat[i][j]<>0 then
        VECT:=ListWithIdenticalEntries(n, 0);
        VECT[i]:=1;
        VECT[j]:=1;
        return VECT;
      fi;
    od;
  od;
end;


DiagonalizeNonDegenerateSymmetricMatrix:=function(SymMat)
  local n, TheBasis, i, BasisOrthogonal, RedInBasis, V, RedMat;
  if SymMat<>TransposedMat(SymMat) then
    Error("The matrix is not symmetric");
  fi;
  n:=Length(SymMat);
  TheBasis:=[];
  for i in [1..n]
  do
    if Length(TheBasis)=0 then
      BasisOrthogonal:=IdentityMat(n);
    else
      BasisOrthogonal:=NullspaceMat(SymMat*TransposedMat(TheBasis));
    fi;
    RedInBasis:=BasisOrthogonal*SymMat*TransposedMat(BasisOrthogonal);
    V:=FindNonIsotropeVector(RedInBasis);
    Add(TheBasis, V*BasisOrthogonal);
  od;
  if Length(TheBasis)=0 then
      RedMat:=[];
  else
      RedMat:=TheBasis*SymMat*TransposedMat(TheBasis);
  fi;
  return rec(RedMat:=RedMat, Transform:=TheBasis);
end;


DiagonalizeSymmetricMatrix:=function(SymMat)
  local NSP, n1, NSP1, RMat1, SymMat2, NSP2, n2, RMat2, i, j, RMat, RedMat, nbZero, nbPlus, nbMinus, strSign;
  n1:=Length(SymMat);
  if SymMat<>TransposedMat(SymMat) then
    Error("The matrix is not symmetric");
  fi;
  NSP1:=NullspaceReduction(SymMat);
  RMat1:=NSP1.Transform;
  SymMat2:=NSP1.NonDegenerate;
  NSP2:=DiagonalizeNonDegenerateSymmetricMatrix(SymMat2);
  n2:=Length(SymMat2);
  RMat2:=NullMat(n1, n1);
  for i in [1..n2]
  do
    for j in [1..n2]
    do
      RMat2[i+n1-n2][j+n1-n2]:=NSP2.Transform[i][j];
    od;
  od;
  for i in [1..n1-n2]
  do
    RMat2[i][i]:=1;
  od;
  RMat:=RMat2*RMat1;
  RedMat:=RMat * SymMat * TransposedMat(RMat);
  nbZero:=0;
  nbPlus:=0;
  nbMinus:=0;
  for i in [1..n1]
  do
    if RedMat[i][i]=0 then
      nbZero:=nbZero+1;
    fi;
    if RedMat[i][i]>0 then
      nbPlus:=nbPlus+1;
    fi;
    if RedMat[i][i]<0 then
      nbMinus:=nbMinus+1;
    fi;
  od;
  strSign:=Concatenation("nbMinus=", String(nbMinus), " nbZero=", String(nbZero), " nbPlus=", String(nbPlus));
  return rec(Transform:=RMat, RedMat:=RedMat, nbZero:=nbZero, nbPlus:=nbPlus, nbMinus:=nbMinus, strSign:=strSign);
end;


GetSetNegativeOrZeroVector:=function(SymMat)
  local TheRec, DiagMat, TheSet, i, n, eVect, fVect, eMax;
  TheRec:=DiagonalizeSymmetricMatrix(SymMat);
  DiagMat:=TheRec.Transform*SymMat*TransposedMat(TheRec.Transform);
  TheSet:=[];
  n:=Length(SymMat);
  for i in [1..n]
  do
    if DiagMat[i][i] <= 0 then
      eVect:=ListWithIdenticalEntries(n,0);
      eVect[i]:=1;
      fVect:=eVect*TheRec.Transform;
      eMax:=Maximum(List(fVect, AbsInt));
      Add(TheSet, fVect/eMax);
    fi;
  od;
  return TheSet;
end;


GetSomeNegativeVector:=function(SymMat)
    local ListNeg, eMult, eVect, eIntVect;
    ListNeg:=GetSetNegativeOrZeroVector(SymMat);
    eMult:=1;
    while(true)
    do
        for eVect in ListNeg
        do
            eIntVect:=List(eMult*eVect, NearestInteger);
            if eIntVect*SymMat*eIntVect < 0 then
                if IndefinitePrint then
                    Print("GetSomeNegativeVector, eMult=", eMult, "\n");
                fi;
                return eIntVect;
            fi;
        od;
        eMult:=eMult+1;
    od;
end;


FindNegativeVect:=function(GramMat)
  return GetSomeNegativeVector(GramMat);
end;
