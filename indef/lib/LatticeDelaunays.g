DelaunayDatabaseManagement:=function()
    local ListDelaunayEXT, ListDelaunayINV, ListDelaunayGroup, ListDelaunayAdjacencies, ListDelaunayStatus, FuncInsertAdjacencies, FuncDelaunayGetNumber, FuncDelaunayGetEXT, FuncDelaunayGetINV, FuncDelaunayGetGroup, FuncDelaunayGetAdjacencies, FuncDelaunayGetStatus, FuncReturnCompleteDescription, FuncInsertDelaunay, FuncDelaunayGetNbEXT, ListDelaunayNbEXT, FuncReturnSingleDelaunayComplete, IsInitialized, GetInitState;
    IsInitialized:=false;
    ListDelaunayEXT:=[];
    ListDelaunayINV:=[];
    ListDelaunayNbEXT:=[];
    ListDelaunayGroup:=[];
    ListDelaunayAdjacencies:=[];
    ListDelaunayStatus:=[];
    FuncInsertDelaunay:=function(TheEXT, TheINV, TheStab)
        local nbDelaunay;
        IsInitialized:=true;
        Add(ListDelaunayINV, TheINV);
        Add(ListDelaunayStatus, "NO");
        Add(ListDelaunayNbEXT, Length(TheEXT));
        Add(ListDelaunayEXT, TheEXT);
        Add(ListDelaunayGroup, TheStab);
        nbDelaunay:=Length(ListDelaunayStatus);
    end;
    FuncInsertAdjacencies:=function(iDel, Adjacencies)
        IsInitialized:=true;
        ListDelaunayAdjacencies[iDel]:=Adjacencies;
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
               FuncDelaunayGetAdjacencies:=FuncDelaunayGetAdjacencies,
               FuncDelaunayGetStatus:=FuncDelaunayGetStatus,
               FuncReturnSingleDelaunayComplete:=FuncReturnSingleDelaunayComplete,
               FuncReturnCompleteDescription:=FuncReturnCompleteDescription);
end;


ComputeDelaunayDecomposition:=function(DataLattice, DataPolyhedral, DelaunayDatabase)
    local n, EXT, FuncInsert, iOrb, IsFinished, EST, Adjacencies, EXTnew, TheStab, BF, ListOrbit, eOrb, iOrbAdj, iOrbSelect, MinSize, nbV, IsFirst, TheAdj, TheTestAdj;
    n:=DataLattice.n;
    BF:=BankRecording(DataPolyhedral.FuncStabilizer, DataPolyhedral.FuncIsomorphy, DataPolyhedral.FuncInvariant, DataPolyhedral.GroupFormalism);
    FuncInsert:=function(EXT)
        local MyInv, iDelaunay, reply, TheStab, TheEXT, TheTest;
        MyInv:=DataLattice.FuncInvariant(DataLattice, EXT);
        if IndefinitePrint then
            Print("|Database|=", DelaunayDatabase.FuncDelaunayGetNumber(), "\n");
        fi;
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
        if IndefinitePrint then
            Print("Find polyhedral object with |EXT|=", Length(EXT), "\n");
        fi;
        TheTest:=DataLattice.KillingDelaunay(EXT, MyInv);
        if TheTest<>false then
            return rec(success:=0, Reason:=TheTest);
        fi;
        TheStab:=DataLattice.FuncStabilizerDelaunay(DataLattice, EXT);
        DelaunayDatabase.FuncInsertDelaunay(EXT, MyInv, TheStab);
        if IndefinitePrint then
            Print("Find Delaunay: ");
            Print(" |V|=", Length(EXT), " ");
            Print(" |LattIsom|=", Order(TheStab.PermutationStabilizer));
            Print("\n");
        fi;
        return rec(success:=1, result:=rec(eBigMat:=IdentityMat(n+1), iDelaunay:=DelaunayDatabase.FuncDelaunayGetNumber()));
    end;
    EXT:=DataLattice.FindDelaunayPolytope();
    EST:=FuncInsert(EXT);
    if EST.success=0 then
        return EST.Reason;
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
            if IndefinitePrint then
                Print("Starting the analysis of Delaunay ", iOrbSelect, " with ", Length(EXT), " vertices\n");
                Print("Beginning the polyhedral computation\n");
            fi;
            ListOrbit:=__ListFacetByAdjacencyDecompositionMethod(EXT, TheStab.PermutationStabilizer, DataPolyhedral, BF);
            if IndefinitePrint then
                Print("   Ending the polyhedral computation, |ListOrbit|=", Length(ListOrbit), "\n");
            fi;
            #
            #
            Adjacencies:=[];
            for iOrbAdj in [1..Length(ListOrbit)]
            do
                if IndefinitePrint then
                    Print("iOrbAdj=", iOrbAdj, "/", Length(ListOrbit), "\n");
                fi;
                eOrb:=ListOrbit[iOrbAdj];
                EXTnew:=DataLattice.FindAdjacentDelaunay(EXT, eOrb);
                TheTestAdj:=DataLattice.KillingAdjacency(EXT, EXTnew);
                if TheTestAdj<>false then
                    return TheTestAdj;
                fi;
                EST:=FuncInsert(EXTnew);
                if EST.success=0 then
                    return EST.Reason;
                fi;
                TheAdj:=EST.result;
                TheAdj.eInc:=eOrb;
                Add(Adjacencies, TheAdj);
            od;
            if IndefinitePrint then
                Print("Adjacency work finished for Orbit ", iOrbSelect, "/", DelaunayDatabase.FuncDelaunayGetNumber(), " orbits\n");
            fi;
            #
            #
            DelaunayDatabase.FuncInsertAdjacencies(iOrbSelect, Adjacencies);
        else
            break;
        fi;
    od;
    if IndefinitePrint then
        Print("Delaunay computation finished\n");
    fi;
    return "all was ok";
end;
