FileDR2:=Filename(DirectoriesPackagePrograms("indefinite"),"dreadnaut");
FileNautyGroupGAP:=Filename(DirectoriesPackagePrograms("indefinite"),"NautyGroupToGAP_sec");
FileNautyIsoOutputGAP:=Filename(DirectoriesPackagePrograms("indefinite"),"NautyIsoOutputToGAP");


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


__PrintPartition:=function(output, ThePartition)
  local nbPart, ThePart, j, i;
  AppendTo(output, "f=[");
  nbPart:=Length(ThePartition);
  for i in [1..nbPart]
  do
    ThePart:=ThePartition[i];
    for j in [1..Length(ThePart)]
    do
      AppendTo(output, ThePart[j]-1);
      if j < Length(ThePart) then
        AppendTo(output, " ");
      fi;
    od;
    if i<nbPart then
      AppendTo(output, "|");
    fi;
  od;
  AppendTo(output, "]\n");
end;


__PrintGraph:=function(output, ListAdjacency)
  local n, i, eV;
  n:=Length(ListAdjacency);
  AppendTo(output, "g\n");
  for i in [1..n]
  do
    AppendTo(output, i-1, " :");
    for eV in ListAdjacency[i]
    do
      AppendTo(output, " ", eV-1);
    od;
    AppendTo(output, ";\n");
  od;
end;


__PrintGraph_Scalable:=function(output, eRecGraph)
  local n, i, eV;
  n:=eRecGraph.n;
  AppendTo(output, "g\n");
  for i in [1..n]
  do
    AppendTo(output, i-1, " :");
    for eV in eRecGraph.GetAdjacent(i)
    do
      AppendTo(output, " ", eV-1);
    od;
    AppendTo(output, ";\n");
  od;
end;


SymmetryGroupVertexColoredGraphAdjList:=function(ListAdjacency, ThePartition)
  local FileNauty, FileDR, FileRead, FileError, n, output, TheGroup;
  FileNauty:=Filename(POLYHEDRAL_tmpdir, "GraphInput");
  FileDR:=Filename(POLYHEDRAL_tmpdir, "GraphDRout");
  FileRead:=Filename(POLYHEDRAL_tmpdir, "GraphRead");
  FileError:=Filename(POLYHEDRAL_tmpdir, "GraphError");
  RemoveFileIfExist(FileNauty);
  RemoveFileIfExist(FileDR);
  RemoveFileIfExist(FileRead);
  RemoveFileIfExist(FileError);
  n:=Length(ListAdjacency);
  output:=OutputTextFile(FileNauty, true);
  AppendTo(output, "n=", n, "\n");
  __PrintPartition(output, ThePartition);
  __PrintGraph(output, ListAdjacency);
  AppendTo(output, "x\n");
  CloseStream(output);
  Exec(FileDR2, " < ", FileNauty, " > ", FileDR, " 2>", FileError);
  if IsExistingFile(FileDR)=false or IsExistingFile(FileError)=false then
      Error("The file FileDR and FileError are missing");
  fi;
  Exec(FileNautyGroupGAP, " < ", FileDR, " > ", FileRead);
  if IsExistingFile(FileRead)=false then
      Error("The file FileRead is missing");
  fi;
  if IsEmptyFile(FileError)=false then
    Error("Nonempty error file in SymmetryGroupColoredGraph");
  fi;
  TheGroup:=ReadAsFunction(FileRead)();
#  Print(NullMat(5));
  RemoveFile(FileNauty);
  RemoveFile(FileDR);
  RemoveFile(FileRead);
  RemoveFile(FileError);
  return TheGroup;
end;


SymmetryGroupVertexColoredGraphAdjList_Scalable:=function(eRecGraph)
  local FileNauty, FileDR, FileRead, FileError, n, output, TheGroup;
  FileNauty:=Filename(POLYHEDRAL_tmpdir, "GraphInput");
  FileDR:=Filename(POLYHEDRAL_tmpdir, "GraphDRout");
  FileRead:=Filename(POLYHEDRAL_tmpdir, "GraphRead");
  FileError:=Filename(POLYHEDRAL_tmpdir, "GraphError");
  output:=OutputTextFile(FileNauty, true);
  AppendTo(output, "n=", eRecGraph.n, "\n");
  __PrintPartition(output, eRecGraph.ThePartition);
  __PrintGraph_Scalable(output, eRecGraph);
  AppendTo(output, "x\n");
  CloseStream(output);
  Exec(FileDR2, " < ", FileNauty, " > ", FileDR, " 2>", FileError);
  Exec(FileNautyGroupGAP, " < ", FileDR, " > ", FileRead);
  if IsEmptyFile(FileError)=false then
    Error("Nonempty error file in SymmetryGroupColoredGraph");
  fi;
  TheGroup:=ReadAsFunction(FileRead)();;
#  Print(NullMat(5));
  RemoveFile(FileNauty);
  RemoveFile(FileDR);
  RemoveFile(FileRead);
  RemoveFile(FileError);
  return TheGroup;
end;


EquivalenceVertexColoredGraphAdjList:=function(ListAdjacency1, ListAdjacency2, ThePartition)
  local FileNauty, FileDR, FileRead, FileError, n, output, TheReply;
  FileNauty:=Filename(POLYHEDRAL_tmpdir, "GraphInput");
  FileDR:=Filename(POLYHEDRAL_tmpdir, "GraphDRout");
  FileRead:=Filename(POLYHEDRAL_tmpdir, "GraphRead");
  FileError:=Filename(POLYHEDRAL_tmpdir, "GraphError");
  n:=Length(ListAdjacency1);
  if n<>Length(ListAdjacency2) then
    return false;
  fi;
  output:=OutputTextFile(FileNauty, true);
  AppendTo(output, "n=", n, "\n");
  __PrintPartition(output, ThePartition);
  __PrintGraph(output, ListAdjacency1);
  AppendTo(output, "c x @\n");
  __PrintGraph(output, ListAdjacency2);
  AppendTo(output, "x ##\n");
  CloseStream(output);
  Exec(FileDR2, " < ", FileNauty, " > ", FileDR, " 2>", FileError);
  Exec(FileNautyIsoOutputGAP, " ", FileDR, " > ", FileRead);
  if IsEmptyFile(FileError)=false then
    Error("Error in EquivalenceVertexColoredGraphAdjList");
  fi;
  TheReply:=ReadAsFunction(FileRead)();
  RemoveFile(FileNauty);
  RemoveFile(FileDR);
  RemoveFile(FileRead);
  RemoveFile(FileError);
  return TheReply;
end;


EquivalenceVertexColoredGraphAdjList_Scalable:=function(eRecGraph1, eRecGraph2)
  local FileNauty, FileDR, FileRead, FileError, n, output, TheReply;
  FileNauty:=Filename(POLYHEDRAL_tmpdir, "GraphInput");
  FileDR:=Filename(POLYHEDRAL_tmpdir, "GraphDRout");
  FileRead:=Filename(POLYHEDRAL_tmpdir, "GraphRead");
  FileError:=Filename(POLYHEDRAL_tmpdir, "GraphError");
  n:=eRecGraph1.n;
  output:=OutputTextFile(FileNauty, true);
  AppendTo(output, "n=", n, "\n");
  __PrintPartition(output, eRecGraph1.ThePartition);
  __PrintGraph_Scalable(output, eRecGraph1);
  AppendTo(output, "c x @\n");
  __PrintGraph_Scalable(output, eRecGraph2);
  AppendTo(output, "x ##\n");
  CloseStream(output);
  Exec(FileDR2, " < ", FileNauty, " > ", FileDR, " 2>", FileError);
  Exec(FileNautyIsoOutputGAP, " ", FileDR, " > ", FileRead);
  if IsEmptyFile(FileError)=false then
    Error("Error in EquivalenceVertexColoredGraphAdjList_Scalable");
  fi;
  TheReply:=ReadAsFunction(FileRead)();
  RemoveFile(FileNauty);
  RemoveFile(FileDR);
  RemoveFile(FileRead);
  RemoveFile(FileError);
  return TheReply;
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


__Method4Partition:=function(korig, n)
  local k;
  k:=__Method4FindTheK(korig);
  return List([1..k], x->[n*(x-1)+1..n*(x-1)+n]);
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


Method4AutomGroupEdgeColoredGraph:=function(DistMat, SetV)
  local TheListAdjacency, ThePartition, korig, n, GRP;
  korig:=Length(SetV);
  n:=Length(DistMat);
  TheListAdjacency:=Method4modelEdgeColoredGraph(DistMat, SetV);
  ThePartition:=__Method4Partition(korig, n);
  GRP:=SymmetryGroupVertexColoredGraphAdjList(TheListAdjacency, ThePartition);
  return SecondReduceGroupAction(GRP, [1..n]);
end;


Method4EquivalenceEdgeColoredGraph:=function(DistMat1, DistMat2, SetV)
  local korig, n, TheListAdjacency1, TheListAdjacency2, ThePartition, TheEquiv;
  korig:=Length(SetV);
  n:=Length(DistMat1);
  TheListAdjacency1:=Method4modelEdgeColoredGraph(DistMat1, SetV);
  TheListAdjacency2:=Method4modelEdgeColoredGraph(DistMat2, SetV);
  ThePartition:=__Method4Partition(korig,n);
  TheEquiv:=EquivalenceVertexColoredGraphAdjList(TheListAdjacency1, TheListAdjacency2, ThePartition);
  if TheEquiv=false then
    return false;
  fi;
  return TheEquiv{[1..n]};
end;


AutomorphismGroupEdgeColoredGraph:=function(DistMat)
  local SetV;
  SetV:=__SetValue(DistMat);
  return Method4AutomGroupEdgeColoredGraph(DistMat, SetV);
end;


IsIsomorphicEdgeColoredGraph:=function(DistMat1, DistMat2)
  local SetV, k, n, Meth2_NBV, Meth3_NBV;
  SetV:=__SetValue(DistMat1);
  Print("SetV1=", __SetValue(DistMat1), " SetV2=", __SetValue(DistMat2), "\n");
  if __SetValue(DistMat2)<>SetV then
    return false;
  fi;
  k:=Length(SetV);
  n:=Length(DistMat1);
  return Method4EquivalenceEdgeColoredGraph(DistMat1, DistMat2, SetV);
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


GetSetColor_Scalable:=function(eRecScalColor)
  local n, SetV, i, eLine, eRecN;
  SetV:=Set([]);
  n:=eRecScalColor.n;
  for i in [1..n]
  do
#    Print("i=", i, "\n");
    eLine:=eRecScalColor.GetLineColor(i);
    SetV:=Union(SetV, Set(eLine));
  od;
  eRecN:=GetFirstSecondVal(SetV);
  return rec(SetV:=SetV, eRecN:=eRecN);
end;


GetFuncListAdjacentMethod4_Scalable:=function(InfC, eRecScalColor)
  local FirstNewVal, SecondNewVal, korig, k, List01Sets, SetVext, n, n2, nExt, GetLineExtended, GetAdjacent, ThePartition;
  FirstNewVal:=InfC.eRecN[1];
  SecondNewVal:=InfC.eRecN[2];
  korig:=Length(InfC.SetV)+2;
  k:=__Method4FindTheK(korig);
  List01Sets:=BuildSet(k, [0,1]);
  SetVext:=Concatenation(InfC.SetV, InfC.eRecN);
  n:=eRecScalColor.n;
  n2:=n+2;
  nExt:=k*n2;
  ThePartition:=__Method4Partition(korig,n2);
  GetLineExtended:=function(i)
    local eLine;
#    Print("k=", k, " korig=", korig, "\n");
#    Print("n=", n, " nExt=", nExt, "\n");
    if i<=n then
      eLine:=eRecScalColor.GetLineColor(i);
      return Concatenation(eLine{[1..i-1]},[0],eLine{[i+1..n]},[eLine[i],FirstNewVal]);
    fi;
    if i=n+1 then
      return Concatenation(List([1..n],x->eRecScalColor.GetScalarColor(x,x)),[0,SecondNewVal]);
    fi;
    if i=n+2 then
      return Concatenation(ListWithIdenticalEntries(n,FirstNewVal),[SecondNewVal,0]);
    fi;
  end;
  GetAdjacent:=function(iExt)
    local iB, iC, i, kPos, ListAdjacent, iK, eAdj, eLine, j, eColor, iColor;
    iB:=iExt-1;
    iC:=iB mod n2;
    i:=iC+1;
    kPos:=1+(iExt-i)/n2;
    ListAdjacent:=[];
    for iK in [1..k]
    do
      if iK<>kPos then
        eAdj:=i+(iK-1)*n2;
        Add(ListAdjacent, eAdj);
      fi;
    od;
    eLine:=GetLineExtended(i);
    for j in [1..n2]
    do
      if j<>i then
        eColor:=eLine[j];
        iColor:=Position(SetVext, eColor);
        if List01Sets[iColor][kPos]=1 then
          eAdj:=j + (kPos-1)*n2;
          Add(ListAdjacent, eAdj);
        fi;
      fi;
    od;
    return ListAdjacent;
  end;
  return rec(GetAdjacent:=GetAdjacent, n:=nExt, ThePartition:=ThePartition);
end;


AutomorphismGroupColoredGraph_Scalable:=function(eRecScalColor)
  local n, InfC, eRecGraph, GRP, NewListGens, eGen, eList, RetGRP;
  n:=eRecScalColor.n;
  InfC:=GetSetColor_Scalable(eRecScalColor);
  eRecGraph:=GetFuncListAdjacentMethod4_Scalable(InfC, eRecScalColor);
  GRP:=SymmetryGroupVertexColoredGraphAdjList_Scalable(eRecGraph);
  NewListGens:=[];
  for eGen in GeneratorsOfGroup(GRP)
  do
    eList:=List([1..n], x->OnPoints(x, eGen));
    Add(NewListGens, PermList(eList));
  od;
  RetGRP:=Group(NewListGens);
  SetSize(RetGRP, Order(GRP));
  return RetGRP;
end;


AutomorphismGroupColoredGraph:=function(ScalarMat)
  local DistMat, NewListGens, GRP, eGen, eList, RetGRP;
  DistMat:=MappedScalarMatrixDistanceMatrix(ScalarMat);
  NewListGens:=[];
  GRP:=AutomorphismGroupEdgeColoredGraph(DistMat);
  for eGen in GeneratorsOfGroup(GRP)
  do
    eList:=List([1..Length(ScalarMat)], x->OnPoints(x, eGen));
    Add(NewListGens, PermList(eList));
  od;
  RetGRP:=Group(NewListGens);
  SetSize(RetGRP, Order(GRP));
  return RetGRP;
end;


IsIsomorphicColoredGraph_Scalable:=function(eRecScalColor1, eRecScalColor2)
  local n, InfC, eRecGraph1, eRecGraph2, test;
  n:=eRecScalColor1.n;
  if n<>eRecScalColor2.n then
    return false;
  fi;
  InfC:=GetSetColor_Scalable(eRecScalColor1);
  if InfC<>GetSetColor_Scalable(eRecScalColor2) then
    return false;
  fi;
  eRecGraph1:=GetFuncListAdjacentMethod4_Scalable(InfC, eRecScalColor1);
  eRecGraph2:=GetFuncListAdjacentMethod4_Scalable(InfC, eRecScalColor2);
  test:=EquivalenceVertexColoredGraphAdjList_Scalable(eRecGraph1, eRecGraph2);
  if test=false then
    return false;
  else
    return test{[1..n]};
  fi;
end;


IsIsomorphicColoredGraph:=function(ScalarMat1, ScalarMat2)
  local DistMat1, DistMat2, test;
  DistMat1:=MappedScalarMatrixDistanceMatrix(ScalarMat1);
  DistMat2:=MappedScalarMatrixDistanceMatrix(ScalarMat2);
  test:=IsIsomorphicEdgeColoredGraph(DistMat1, DistMat2);
  if test=false then
    return false;
  else
    return test{[1..Length(ScalarMat1)]};
  fi;
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
