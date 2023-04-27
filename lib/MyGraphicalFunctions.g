__SetValue:=function(DistMat)
  local ListSet, n, i, j;
  ListSet:=[];
  n:=Length(DistMat);
  for i in [1..n-1]
  do
    Append(ListSet, DistMat[i]{[i+1..n]});
  od;
  return Set(ListSet);
end;


__SetValue_ScalarMat:=function(ScalarMat)
  local ListSet, n, i, j;
  ListSet:=Set([]);
  n:=Length(ScalarMat);
  for i in [1..n]
  do
    ListSet:=Union(ListSet, Set(ScalarMat[i]{[i..n]}));
  od;
  return ListSet;
end;


__Method4FindTheK:=function(korig)
  local k;
  k:=1;
  while(true)
  do
    if korig < 2^k then
      break;
    fi;
    k:=k+1;
  od;
  return k;
end;


# this procedure Build the Set:  Seto x Seto x .... x Seto
BuildSet:=function(n, Seto)
  local DO, i, iCol, U, V,C, eVal;
  DO:=[[]];
  for iCol in [1..n]
  do
    C:=ShallowCopy(DO);
    DO:=ShallowCopy([]);
    for i in [1..Length(C)]
    do
      for eVal in Seto
      do
        U:=ShallowCopy(C[i]);
        Add(U, eVal);
        Add(DO, U);
      od;
    od;
  od;
  return DO;
end;


Method4modelEdgeColoredGraph:=function(DistMat, SetV)
  local korig, k, List01Sets, n, ListAdjacency, iVert, i, j, iColor, eVal, jVert, iCol;
  korig:=Length(SetV);
  k:=__Method4FindTheK(korig);
  List01Sets:=BuildSet(k, [0,1]);
  n:=Length(DistMat);
  ListAdjacency:=[];
  for i in [1..k*n]
  do
    Add(ListAdjacency, []);
  od;
  for iVert in [1..n]
  do
    for i in [1..k-1]
    do
      for j in [i+1..k]
      do
        Add(ListAdjacency[n*(i-1)+iVert], n*(j-1)+iVert);
        Add(ListAdjacency[n*(j-1)+iVert], n*(i-1)+iVert);
      od;
    od;
  od;
  for iColor in [1..Length(SetV)]
  do
    eVal:=SetV[iColor];
    for iVert in [1..n-1]
    do
      for jVert in [iVert+1..n]
      do
        if DistMat[iVert][jVert]=eVal then
          for iCol in [1..k]
          do
            if List01Sets[iColor][iCol]=1 then
              Add(ListAdjacency[n*(iCol-1)+iVert], n*(iCol-1)+jVert);
              Add(ListAdjacency[n*(iCol-1)+jVert], n*(iCol-1)+iVert);
            fi;
          od;
        fi;
      od;
    od;
  od;
  return ListAdjacency;
end;


GetFirstSecondVal:=function(SetV)
  local FirstNewVal, SecondNewVal;
  FirstNewVal:=0;
  while(true)
  do
    if Position(SetV, FirstNewVal)=fail then
      break;
    fi;
    FirstNewVal:=FirstNewVal+1;
  od;
  SecondNewVal:=FirstNewVal+1;
  while(true)
  do
    if Position(SetV, SecondNewVal)=fail then
      break;
    fi;
    SecondNewVal:=SecondNewVal+1;
  od;
  return [FirstNewVal, SecondNewVal];
end;


#
# we work in the framework of scalar products
# while the dreadnaut formalism uses distance matrices
# i.e. with 0 on the diagonal. We ought to take care of that
MappedScalarMatrixDistanceMatrix:=function(ScalarMat)
  local DistMat, SetV, eRecN, n, n2, i, j, FirstNewVal, SecondNewVal;
  n:=Length(ScalarMat);
  n2:=n+2;
  SetV:=__SetValue_ScalarMat(ScalarMat);
  eRecN:=GetFirstSecondVal(SetV);
  FirstNewVal:=eRecN[1];
  SecondNewVal:=eRecN[2];
  DistMat:=NullMat(n2, n2);
  for i in [1..n-1]
  do
    for j in [i+1..n]
    do
      DistMat[i][j]:=ScalarMat[i][j];
      DistMat[j][i]:=ScalarMat[i][j];
    od;
  od;
  for i in [1..n]
  do
    DistMat[i][n+1]:=ScalarMat[i][i];
    DistMat[n+1][i]:=ScalarMat[i][i];
    DistMat[i][n+2]:=FirstNewVal;
    DistMat[n+2][i]:=FirstNewVal;
  od;
  DistMat[n+1][n+2]:=SecondNewVal;
  DistMat[n+2][n+1]:=SecondNewVal;
  return DistMat;
end;


FindShortestPath:=function(Gra, x, y)
  local FuncNextVert, Dist, ePath, CurrentVert;
  FuncNextVert:=function(VVert, d)
    local u;
    for u in Adjacency(Gra, VVert)
    do
      if Distance(Gra, u, y)=d-1 then
        return u;
      fi;
    od;
  end;
  Dist:=Distance(Gra, x, y);
  if Dist=-1 then
    Error("Cannot find shortest path when vertices are not connected");
  fi;
  ePath:=[x];
  CurrentVert:=x;
  while(Dist>0)
  do
    CurrentVert:=FuncNextVert(CurrentVert, Dist);
    Add(ePath, CurrentVert);
    Dist:=Dist-1;
  od;
  return ePath;
end;
