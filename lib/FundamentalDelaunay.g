

FindDelaunayPolytope_Rational:=function(GramMat)
  local GramMatInt, n, ListCosetRed, ListCosetDiff, ListRelevantPoints, i, V, DefiningInequality, TheRandomDirection, TheLP, eVect, TheNorm, TheCVP, ListInequalities, eEnt, RetEXT;
  Print("Beginning of FindDelaunayPolytope\n");
  GramMatInt:=RemoveFractionMatrix(GramMat);
  n:=Length(GramMat);
  ListRelevantPoints:=[];
  for i in [1..n]
  do
    V:=ListWithIdenticalEntries(n,0);
    V[i]:=1;
    Add(ListRelevantPoints, ShallowCopy(V));
    V[i]:=-1;
    Add(ListRelevantPoints, ShallowCopy(V));
  od;
  DefiningInequality:=function(eVect)
    return Concatenation([eVect*GramMatInt*eVect/2], -eVect*GramMatInt);
  end;
  TheRandomDirection:=FuncRandomDirection(n, 10);
  while(true)
  do
    ListInequalities:=List(ListRelevantPoints, DefiningInequality);
    TheLP:=LinearProgramming_Rational(ListInequalities, TheRandomDirection);
    eVect:=ListWithIdenticalEntries(n,0);
    for eEnt in TheLP.primal_solution
    do
      eVect[eEnt[1]]:=eEnt[2];
    od;
    TheNorm:=eVect*GramMatInt*eVect;
    Print("TheNorm=", TheNorm, "\n");
    TheCVP:=CVPVallentinProgram_Rational(GramMatInt, eVect);
    if TheCVP.TheNorm=TheNorm then
      RetEXT:=List(TheCVP.ListVect, x->Concatenation([1], x));
      Print("Ending of FindDelaunayPolytope |RetEXT|=", Length(RetEXT), "\n");
      return RetEXT;
    fi;
    if TheNorm < TheCVP.TheNorm then
      Error("Inconsistency in norm computation");
    fi;
    for eVect in TheCVP.ListVect
    do
      Add(ListRelevantPoints, eVect);
    od;
    Print("|ListRelevantPoints|=", Length(ListRelevantPoints), "\n");
  od;
end;


FindDelaunayPolytope_QN:=function(Nval, GramMat)
  local n, ListCosetRed, ListCosetDiff, ListRelevantPoints, i, V, DefiningInequality, TheRandomDirection, TheLP, eVect, TheNorm, TheCVP, ListInequalities, eEnt, RetEXT;
  Print("Beginning of FindDelaunayPolytope_QN\n");
  n:=Length(GramMat);
  ListRelevantPoints:=[];
  for i in [1..n]
  do
    V:=ListWithIdenticalEntries(n,0);
    V[i]:=1;
    Add(ListRelevantPoints, ShallowCopy(V));
    V[i]:=-1;
    Add(ListRelevantPoints, ShallowCopy(V));
  od;
  DefiningInequality:=function(eVect)
    return Concatenation([eVect*GramMat*eVect/2], -eVect*GramMat);
  end;
  TheRandomDirection:=FuncRandomDirection(n, 10);
  while(true)
  do
    ListInequalities:=List(ListRelevantPoints, DefiningInequality);
    TheLP:=LinearProgramming_QN(Nval, ListInequalities, TheRandomDirection);
    eVect:=ListWithIdenticalEntries(n,0);
    for eEnt in TheLP.primal_solution
    do
      eVect[eEnt[1]]:=eEnt[2];
    od;
    TheNorm:=eVect*GramMat*eVect;
#    Print("TheNorm=", TheNorm, "\n");
    TheCVP:=CVPVallentinProgram_QN(Nval, GramMat, eVect);
    if TheCVP.TheNorm=TheNorm then
      RetEXT:=List(TheCVP.ListVect, x->Concatenation([1], x));
      Print("Ending of FindDelaunayPolytope_QN |RetEXT|=", Length(RetEXT), "\n");
      return RetEXT;
    fi;
    if QN_IsPositive(Nval, TheNorm-TheCVP.TheNorm)=false then
      Error("Inconsistency in norm computation");
    fi;
    for eVect in TheCVP.ListVect
    do
      Add(ListRelevantPoints, eVect);
    od;
    Print("|ListRelevantPoints|=", Length(ListRelevantPoints), "\n");
  od;
end;


FindDelaunayPolytope:=function(GramMat)
  local Nval;
  if IsMatrixRational(GramMat)=true then
    return FindDelaunayPolytope_Rational(GramMat);
  fi;
  for Nval in [2,5]
  do
    if QN_IsMatrix(Nval, GramMat)=true then
      return FindDelaunayPolytope_QN(Nval, GramMat);
    fi;
  od;
  Error("You have to build your own arithmetic");
end;
