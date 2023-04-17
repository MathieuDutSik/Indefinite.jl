FindDelaunayPolytope:=function(GramMat)
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
    TheLP:=LinearProgramming(ListInequalities, TheRandomDirection);
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
