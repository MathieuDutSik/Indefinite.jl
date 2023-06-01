# this is the program for Bank management
# it does not assume a fixed EXT set, so as to deal with
# polyhedral subdivision.
# FuncStabilizer take as argument EXT and ListInc and return a stabilizer
# of EXT. N
# FuncIsomorphy(EXT1,  EXT2) and return a
# false if there is no isomorphy and a permutation ePerm where the ith
# element is the position in ListInc2 of ListInc1[i].
#
# the system is not designed for a large number of accounts, just for
# a few hundreds accounts.
BankRecording:=function(FuncStabilizer, FuncIsomorphy, FuncInvariant, GroupFormalism)
  local FuncRetrieveObject, FuncCreateAccount, FuncStab, WRL, ListCompleteInformation, ListInvariant, eInv, MinNbVert, nbRec, iRec, eInfoEXT, nbVert;
  nbRec:=0;
  #
  MinNbVert:=-1;
  WRL:=[];
  #
  ListCompleteInformation:=[];
  ListInvariant:=[];
  FuncRetrieveObject:=function(EXT, GivenSymmetry)
    local eChar, iAccount, eTransform, CompleteAccount, EXTaccount, ListObjectAccount, TransListObject, GRPaccount, TheGrp, NewGrp, eInc, RPL, eInvariant;
    if Length(EXT) < MinNbVert then
      return false;
    fi;
    eChar:=GroupFormalism.BankKeyInformation(EXT, GivenSymmetry);
    eInvariant:=FuncInvariant(EXT);
    for iAccount in [1..nbRec]
    do
      if eInvariant=ListInvariant[iAccount] then
        eInfoEXT:=WRL[iAccount];
        eTransform:=FuncIsomorphy(GroupFormalism.BankGetForIsom(eInfoEXT), GroupFormalism.BankGetForIsom(eChar));
        if eTransform<>false then
          CompleteAccount:=ListCompleteInformation[iAccount];
          EXTaccount:=GroupFormalism.BankGetVertexSet(eInfoEXT, CompleteAccount);
          ListObjectAccount:=GroupFormalism.BankGetListObject(CompleteAccount);
          TransListObject:=GroupFormalism.TransformIncidenceList(EXTaccount, EXT, eTransform, ListObjectAccount);
          GRPaccount:=GroupFormalism.BankGetGroup(eInfoEXT, CompleteAccount);
          TheGrp:=GroupFormalism.GroupConjugacy(GRPaccount, eTransform);
          if GroupFormalism.IsSubgroup(TheGrp, GivenSymmetry)=true then
            return rec(GRP:=TheGrp, ListOrbitFacet:=TransListObject);
          else
            NewGrp:=GroupFormalism.GroupUnion(TheGrp, GivenSymmetry);
            RPL:=GroupFormalism.OrbitGroupFormalism(EXT, NewGrp);
            for eInc in TransListObject
            do
              RPL.FuncInsert(eInc);
            od;
            return rec(GRP:=NewGrp, ListOrbitFacet:=RPL.FuncListOrbitIncidence());
          fi;
        fi;
      fi;
    od;
    return false;
  end;
  FuncCreateAccount:=function(EXT, GroupExt, ListObject)
    local CompleteInfo, eInvariant, eInfoEXT, nbVert;
    eInfoEXT:=GroupFormalism.BankKeyInformation(EXT, GroupExt);
    nbVert:=Length(eInfoEXT.EXT);
    if nbRec=0 then
      MinNbVert:=nbVert;
    else
      MinNbVert:=Minimum(MinNbVert, nbVert);
    fi;
    eInvariant:=FuncInvariant(EXT);
    nbRec:=nbRec+1;
    #
    CompleteInfo:=GroupFormalism.BankCompleteInformation(EXT, GroupExt, ListObject);
    Add(ListCompleteInformation, CompleteInfo);
    Add(WRL, eInfoEXT);
    Add(ListInvariant, eInvariant);
  end;
  return rec(FuncStabilizer:=FuncStabilizer, FuncCreateAccount:=FuncCreateAccount, FuncRetrieveObject:=FuncRetrieveObject);
end;


__ListFacetByAdjacencyDecompositionMethod:=function(EXT, GivenSymmetry, Data, BankFormalism)
  local RECListOrbit, IsFinished, eOrb, eInc, Ladj, SelectedOrbit, jOrb, MiniINCD, RPL, RPLift, testBank, OrdGRP, TheDim, WorkingSymGroup, LRES, NewData, RedStab, TestNeedMoreSymmetry, ReturnedList, eSetUndone, nbUndone, testSym, fInc;
  testBank:=BankFormalism.FuncRetrieveObject(EXT, GivenSymmetry);
  if testBank<>false then
    return Data.GroupFormalism.LiftingOrbits(EXT, testBank.ListOrbitFacet, GivenSymmetry, testBank.GRP);
  fi;
  # we would like to use IsBankSave but it is not possible with EllaspedTime
  TestNeedMoreSymmetry:=function(EXT)
    if Length(EXT)> RankMat(EXT)+4 then
      return true;
    else
      return false;
    fi;
  end;
  if IsBound(Data.TestNeedMoreSymmetry)=true then
    testSym:=Data.TestNeedMoreSymmetry(EXT);
  else
    testSym:=TestNeedMoreSymmetry(EXT);
  fi;
  if testSym=true then
    WorkingSymGroup:=Data.GroupFormalism.GroupUnion(BankFormalism.FuncStabilizer(EXT), GivenSymmetry);
  else
    WorkingSymGroup:=GivenSymmetry;
  fi;
  OrdGRP:=Data.GroupFormalism.Order(WorkingSymGroup);
  if Data.IsRespawn(OrdGRP, EXT, Data.TheDepth)=false then
    ReturnedList:=Data.DualDescriptionFunction(EXT, Data.GroupFormalism.ToPermGroup(EXT, GivenSymmetry));
    if Data.IsBankSave(OrdGRP, EXT, Data.TheDepth)=true then
#      Print("Before FuncCreateAccount\n");
      BankFormalism.FuncCreateAccount(EXT, GivenSymmetry, ReturnedList);
#      Print("After FuncCreateAccount\n");
    fi;
  else
    TheDim:=RankMat(EXT)-1;
#    Print("RESPAWN a new ADM computation |GRP|=", OrdGRP, " TheDim=", TheDim, " |EXT|=", Length(EXT), "\n");
    RPL:=Data.GroupFormalism.OrbitGroupFormalism(EXT, WorkingSymGroup);
#    Print("nbOrbit=", RPL.FuncNumberOrbit(), "\n");
    if RPL.FuncNumberOrbit()=0 then
      for eInc in Data.GetInitialRays(EXT, 10)
      do
        RPL.FuncInsert(eInc);
      od;
    fi;
#    Print("Running the adjacency method recursively\n");
    while(true)
    do
      eSetUndone:=RPL.ComputeIntersectionUndone();
      nbUndone:=RPL.FuncNumberUndone();
      if RPL.FuncNumberOrbitDone()>0 then
        if nbUndone<=TheDim-1 or Length(eSetUndone)>0 then
#          Print("End of computation, nbObj=", RPL.FuncNumber(), " nbOrbit=", RPL.FuncNumberOrbit(), " nbUndone=", nbUndone, " |eSetUndone|=", Length(eSetUndone), " Depth=", Data.TheDepth, " |EXT|=", Length(EXT), "\n");
          break;
        fi;
      fi;
      SelectedOrbit:=RPL.FuncGetMinimalUndoneOrbit();
      eInc:=RPL.FuncRecord(SelectedOrbit);
#      Print("\n");
      RedStab:=Data.GroupFormalism.Stabilizer(EXT, WorkingSymGroup, eInc);
#      Print("Considering orbit ", SelectedOrbit, " |inc|=", Length(eInc), " Depth=", Data.TheDepth, " |stab|=", Order(RedStab), " dim=", TheDim, "\n");
      RPLift:=__ProjectionLiftingFramework(EXT, eInc);
      NewData:=ShallowCopy(Data);
      NewData.TheDepth:=NewData.TheDepth+1;
      Ladj:=__ListFacetByAdjacencyDecompositionMethod(EXT{eInc}, RedStab, NewData, BankFormalism);
#      Print("We treat ", Length(Ladj), " orbits\n");
      for fInc in Ladj
      do
        RPL.FuncInsert(RPLift.FuncLift(fInc));
      od;
      RPL.FuncPutOrbitAsDone(SelectedOrbit);
      nbUndone:=RPL.FuncNumberUndone();
#      Print("We have ", RPL.FuncNumberOrbit(), " orbits");
#      Print("  Nr treated=", RPL.FuncNumberOrbitDone(), " orbits");
#      Print("  nbUndone=", nbUndone);
#      Print("\n");
    od;
    LRES:=RPL.FuncListOrbitIncidence();
    ReturnedList:=Data.GroupFormalism.LiftingOrbits(EXT, LRES, GivenSymmetry, WorkingSymGroup);
#    Print("LRES=", Length(LRES), " |ReturnedList|=", Length(ReturnedList), " |GivenSymmetry|=", Order(GivenSymmetry), " |WorkingSymGroup|=", Order(WorkingSymGroup), "\n");
    if Data.IsBankSave(OrdGRP, EXT, Data.TheDepth)=true then
#      Print("Before FuncCreateAccount\n");
      BankFormalism.FuncCreateAccount(EXT, WorkingSymGroup, LRES);
#      Print("After FuncCreateAccount\n");
    fi;
  fi;
  return ReturnedList;
end;
