DelaunayDatabaseManagement:=function()
  local ListDelaunayEXT, ListDelaunayINV, ListDelaunayGroup, ListDelaunayAdjacencies, ListDelaunayStatus, FuncInsertAdjacencies, FuncDelaunayGetNumber, FuncDelaunayGetEXT, FuncDelaunayGetINV, FuncDelaunayGetGroup, FuncDelaunayGetAdjacencies, FuncDelaunayGetStatus, FuncReturnCompleteDescription, FuncInsertDelaunay, FuncDestroyDatabase, FuncDelaunayGetNbEXT, ListDelaunayNbEXT, FuncReturnSingleDelaunayComplete, IsInitialized, GetInitState, PathDelaunayPOLY;
  IsInitialized:=false;
  ListDelaunayEXT:=[];
  ListDelaunayINV:=[];
  ListDelaunayNbEXT:=[];
  ListDelaunayGroup:=[];
  ListDelaunayAdjacencies:=[];
  ListDelaunayStatus:=[];
  FuncInsertDelaunay:=function(TheEXT, TheINV, TheStab)
    local nbDelaunay, FileDelaunayEXT, FileDelaunayINV, FileGroup;
    IsInitialized:=true;
    Add(ListDelaunayINV, TheINV);
    Add(ListDelaunayStatus, "NO");
    Add(ListDelaunayNbEXT, Length(TheEXT));
    if MemorySave=false then
      Add(ListDelaunayEXT, TheEXT);
      Add(ListDelaunayGroup, TheStab);
    fi;
    nbDelaunay:=Length(ListDelaunayStatus);
  end;
  FuncInsertAdjacencies:=function(iDel, Adjacencies)
    local FileDelaunayAdd;
    IsInitialized:=true;
    if MemorySave=false then
      ListDelaunayAdjacencies[iDel]:=Adjacencies;
    fi;
    ListDelaunayStatus[iDel]:="YES";
  end;
  FuncDelaunayGetNumber:=function()
    return Length(ListDelaunayINV);
  end;
  FuncDelaunayGetNbEXT:=function(iOrb)
    return ListDelaunayNbEXT[iOrb];
  end;
  FuncDelaunayGetEXT:=function(iOrb)
    return ListDelaunayEXT[iOrb];
  end;
  FuncDelaunayGetINV:=function(iOrb)
    return ListDelaunayINV[iOrb];
  end;
  FuncDelaunayGetGroup:=function(iOrb)
    return ListDelaunayGroup[iOrb];
  end;
  FuncDelaunayGetStatus:=function(iOrb)
    return ListDelaunayStatus[iOrb];
  end;
  FuncDelaunayGetAdjacencies:=function(iOrb)
    return ListDelaunayAdjacencies[iOrb];
  end;
  FuncReturnSingleDelaunayComplete:=function(iOrb)
    return rec(EXT:=FuncDelaunayGetEXT(iOrb),
               Linv:=FuncDelaunayGetINV(iOrb),
               TheStab:=FuncDelaunayGetGroup(iOrb),
               Adjacencies:=FuncDelaunayGetAdjacencies(iOrb));
  end;
  FuncReturnCompleteDescription:=function()
    local ListOrbitDelaunay, iOrb;
    ListOrbitDelaunay:=[];
    for iOrb in [1..FuncDelaunayGetNumber()]
    do
      Add(ListOrbitDelaunay, FuncReturnSingleDelaunayComplete(iOrb));
    od;
    return ListOrbitDelaunay;
  end;
  GetInitState:=function()
    return IsInitialized;
  end;
  return rec(GetInitState:=GetInitState,
             FuncInsertDelaunay:=FuncInsertDelaunay,
             FuncInsertAdjacencies:=FuncInsertAdjacencies,
             FuncDelaunayGetNbEXT:=FuncDelaunayGetNbEXT,
             FuncDelaunayGetNumber:=FuncDelaunayGetNumber,
             FuncDelaunayGetEXT:=FuncDelaunayGetEXT,
             FuncDelaunayGetINV:=FuncDelaunayGetINV,
             FuncDelaunayGetGroup:=FuncDelaunayGetGroup,
             FuncDestroyDatabase:=FuncDestroyDatabase,
             FuncDelaunayGetAdjacencies:=FuncDelaunayGetAdjacencies,
             FuncDelaunayGetStatus:=FuncDelaunayGetStatus,
             FuncReturnSingleDelaunayComplete:=FuncReturnSingleDelaunayComplete,
             FuncReturnCompleteDescription:=FuncReturnCompleteDescription);
end;


ComputeDelaunayDecomposition:=function(DataLattice, DataPolyhedral, DelaunayDatabase)
  local n, EXT, FuncInsert, iOrb, IsFinished, EST, Adjacencies, EXTnew, TheStab, BF, FilePolyhedralOrb, ListOrbit, eOrb, iOrbAdj, iOrbSelect, ThePath, FileSingleAdjacency, BankPath, FuncClearComputation, MinSize, nbV, IsFirst, TheAdj, TheTestAdj;
  ThePath:=DataLattice.PathPermanent;
  BankPath:=Concatenation(ThePath, "TheBank/");
  if DataLattice.Saving=true then
    Exec("mkdir -p ", BankPath);
  fi;
  n:=DataLattice.n;
  FuncClearComputation:=function()
    BF.FuncClearAccount();
    DelaunayDatabase.FuncDestroyDatabase();
  end;
  BF:=BankRecording(rec(Saving:=DataLattice.Saving, BankPath:=BankPath), DataPolyhedral.FuncStabilizer, DataPolyhedral.FuncIsomorphy, DataPolyhedral.FuncInvariant, DataPolyhedral.GroupFormalism);
  FuncInsert:=function(EXT)
    local MyInv, iDelaunay, reply, TheStab, TheEXT, TheTest;
    MyInv:=DataLattice.FuncInvariant(DataLattice, EXT);
    Print("|Database|=", DelaunayDatabase.FuncDelaunayGetNumber(), "\n");
    for iDelaunay in [1..DelaunayDatabase.FuncDelaunayGetNumber()]
    do
      if MyInv=DelaunayDatabase.FuncDelaunayGetINV(iDelaunay) then
        TheEXT:=DelaunayDatabase.FuncDelaunayGetEXT(iDelaunay);
        TheStab:=DelaunayDatabase.FuncDelaunayGetGroup(iDelaunay);
        reply:=DataLattice.FuncIsomorphismDelaunay(DataLattice, TheEXT, EXT, TheStab);
        if reply<>false then
          return rec(success:=1, result:=rec(eBigMat:=reply, iDelaunay:=iDelaunay));
        fi;
      fi;
    od;
    Print("Find polyhedral object with |EXT|=", Length(EXT), "\n");
    TheTest:=DataLattice.KillingDelaunay(EXT, MyInv);
    if TheTest<>false then
      return rec(success:=0, Reason:=TheTest);
    fi;
    TheStab:=DataLattice.FuncStabilizerDelaunay(DataLattice, EXT);
    DelaunayDatabase.FuncInsertDelaunay(EXT, MyInv, TheStab);
    Print("Find Delaunay: ");
    Print(" |V|=", Length(EXT), " ");
    Print(" |LattIsom|=", Order(TheStab.PermutationStabilizer));
    Print("\n");
    return rec(success:=1, result:=rec(eBigMat:=IdentityMat(n+1), iDelaunay:=DelaunayDatabase.FuncDelaunayGetNumber()));
  end;
  DelaunayDatabase.Recover();
  if DelaunayDatabase.FuncDelaunayGetNumber()=0 then
    EXT:=DataLattice.FindDelaunayPolytope();
    EST:=FuncInsert(EXT);
    if EST.success=0 then
      FuncClearComputation();
      return EST.Reason;
    fi;
  fi;
  while(true)
  do
    IsFinished:=true;
    IsFirst:=true;
    MinSize:=0;
    iOrbSelect:=-1;
    for iOrb in [1..DelaunayDatabase.FuncDelaunayGetNumber()]
    do
      if DelaunayDatabase.FuncDelaunayGetStatus(iOrb)="NO" then
        nbV:=DelaunayDatabase.FuncDelaunayGetNbEXT(iOrb);
        if nbV < MinSize or IsFirst=true then
          MinSize:=nbV;
          iOrbSelect:=iOrb;
          IsFinished:=false;
        fi;
        IsFirst:=false;
      fi;
    od;
    if IsFinished=false then
      TheStab:=DelaunayDatabase.FuncDelaunayGetGroup(iOrbSelect);
      EXT:=DelaunayDatabase.FuncDelaunayGetEXT(iOrbSelect);
      Print("Starting the analysis of Delaunay ", iOrbSelect, " with ", Length(EXT), " vertices\n");
      Print("Beginning the polyhedral computation\n");
      FilePolyhedralOrb:=Concatenation(ThePath, "ListPOLY/PolyhedralSave", String(iOrbSelect));
      ListOrbit:=ComputeAndSavePlusTouchPlusTest(FilePolyhedralOrb, x->__ListFacetByAdjacencyDecompositionMethod(EXT, TheStab.PermutationStabilizer, DataPolyhedral, BF), DataLattice.Saving);
      Print("   Ending the polyhedral computation, |ListOrbit|=", Length(ListOrbit), "\n");
      #
      #
      Adjacencies:=[];
      for iOrbAdj in [1..Length(ListOrbit)]
      do
          Print("iOrbAdj=", iOrbAdj, "/", Length(ListOrbit), "\n");
          eOrb:=ListOrbit[iOrbAdj];
          EXTnew:=DataLattice.FindAdjacentDelaunay(EXT, eOrb);
          TheTestAdj:=DataLattice.KillingAdjacency(EXT, EXTnew);
          if TheTestAdj<>false then
            FuncClearComputation();
            return TheTestAdj;
          fi;
          EST:=FuncInsert(EXTnew);
          if EST.success=0 then
            FuncClearComputation();
            return EST.Reason;
          fi;
          TheAdj:=EST.result;
          TheAdj.eInc:=eOrb;
          Add(Adjacencies, TheAdj);
      od;
      Print("Adjacency work finished for Orbit ", iOrbSelect, "/", DelaunayDatabase.FuncDelaunayGetNumber(), " orbits\n");
      #
      #
      DelaunayDatabase.FuncInsertAdjacencies(iOrbSelect, Adjacencies);
    else
      break;
    fi;
  od;
  Print("Delaunay computation finished\n");
  return "all was ok";
end;
