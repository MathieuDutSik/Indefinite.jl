

Bank_GetSize:=function(ThePrefix)
  local nbRec, eFileEXT, eFileFAC;
  nbRec:=0;
  while(true)
  do
    eFileEXT:=Concatenation(ThePrefix, "AccountEXT_", String(nbRec+1));
    eFileFAC:=Concatenation(ThePrefix, "AccountFAC_", String(nbRec+1));
    if IsExistingFilePlusTouch(eFileEXT)=false or IsExistingFilePlusTouch(eFileFAC)=false then
      break;
    fi;
    nbRec:=nbRec+1;
  od;
  return nbRec;
end;


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
BankRecording:=function(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, GroupFormalism)
  local FuncRetrieveObject, FuncCreateAccount, FuncClearAccount, FuncStab, FileBankRecord, WRL, ListCompleteInformation, ListInvariant, eInv, MinNbVert, nbRec, iRec, FileNameEXT, FileNameINV, eInfoEXT, nbVert;
  if IsCorrectPath(DataBank.BankPath)=false then
    Print("Variable DataBank.BankPath=", DataBank.BankPath, " is incorrect\n");
    Error("It should finish with a /");
  fi;
  if IsDirectoryPath(DataBank.BankPath)=false and DataBank.Saving=true then
    Error("Directory DataBank.BankPath=", DataBank.BankPath, " is nonexistent");
  fi;
  if DataBank.Saving then
    nbRec:=Bank_GetSize(DataBank.BankPath);
  else
    nbRec:=0;
  fi;
#  Print("Setting up database. Right now nbRec=", nbRec, "\n");
  #
  MinNbVert:=-1;
  if DataBank.Saving then
    for iRec in [1..nbRec]
    do
      FileNameEXT:=Concatenation(DataBank.BankPath, "AccountEXT_", String(iRec));
      eInfoEXT:=ReadAsFunction(FileNameEXT)();
      nbVert:=Length(GroupFormalism.BankGetForIsom(eInfoEXT));
      if iRec=1 then
        MinNbVert:=nbVert;
      else
        MinNbVert:=Minimum(MinNbVert, nbVert);
      fi;
    od;
  else
    WRL:=[];
  fi;
#  Print("At database loading. MinNbVert=", MinNbVert, "\n");
  #
  if DataBank.Saving=false then
    ListCompleteInformation:=[];
  fi;
  #
  ListInvariant:=[];
  for iRec in [1..nbRec]
  do
    FileNameINV:=Concatenation(DataBank.BankPath, "AccountINV_", String(iRec));
    if IsExistingFilePlusTouch(FileNameINV) then
      eInv:=ReadAsFunction(FileNameINV)();
    else
      FileNameEXT:=Concatenation(DataBank.BankPath, "AccountEXT_", String(iRec));
      eInfoEXT:=ReadAsFunction(FileNameEXT)();
      eInv:=FuncInvariant(GroupFormalism.BankGetForIsom(eInfoEXT));
      SaveDataToFilePlusTouch(FileNameINV, eInv);
    fi;
    Add(ListInvariant, eInv);
  od;
#  Print("Invariants have been loaded\n");
  FuncRetrieveObject:=function(EXT, GivenSymmetry)
    local eChar, iAccount, eTransform, FileName, CompleteAccount, EXTaccount, ListObjectAccount, TransListObject, GRPaccount, TheGrp, NewGrp, eInc, RPL, eInvariant;
    if Length(EXT) < MinNbVert then
      return false;
    fi;
    eChar:=GroupFormalism.BankKeyInformation(EXT, GivenSymmetry);
    eInvariant:=FuncInvariant(EXT);
#    Print("Beginning search in |database|=", nbRec, "\n");
    for iAccount in [1..nbRec]
    do
      if eInvariant=ListInvariant[iAccount] then
        if DataBank.Saving=true then
          FileNameEXT:=Concatenation(DataBank.BankPath, "AccountEXT_", String(iAccount));
          eInfoEXT:=ReadAsFunction(FileNameEXT)();
        else
          eInfoEXT:=WRL[iAccount];
        fi;
        eTransform:=FuncIsomorphy(GroupFormalism.BankGetForIsom(eInfoEXT), GroupFormalism.BankGetForIsom(eChar));
        if eTransform<>false then
          if DataBank.Saving=true then
            FileName:=Concatenation(DataBank.BankPath, "AccountFAC_", String(iAccount));
            CompleteAccount:=ReadAsFunction(FileName)();
          else
            CompleteAccount:=ListCompleteInformation[iAccount];
          fi;
          EXTaccount:=GroupFormalism.BankGetVertexSet(eInfoEXT, CompleteAccount);
          ListObjectAccount:=GroupFormalism.BankGetListObject(CompleteAccount);
          TransListObject:=GroupFormalism.TransformIncidenceList(EXTaccount, EXT, eTransform, ListObjectAccount);
          GRPaccount:=GroupFormalism.BankGetGroup(eInfoEXT, CompleteAccount);
          TheGrp:=GroupFormalism.GroupConjugacy(GRPaccount, eTransform);
          if GroupFormalism.IsSubgroup(TheGrp, GivenSymmetry)=true then
            return rec(GRP:=TheGrp, ListOrbitFacet:=TransListObject);
          else
            NewGrp:=GroupFormalism.GroupUnion(TheGrp, GivenSymmetry);
            RPL:=GroupFormalism.OrbitGroupFormalism(EXT, NewGrp, "/irrelevant/", false);
            for eInc in TransListObject
            do
              RPL.FuncInsert(eInc);
            od;
            return rec(GRP:=NewGrp, ListOrbitFacet:=RPL.FuncListOrbitIncidence());
          fi;
        fi;
      fi;
    od;
#    Print("Ending search in database\n");
    return false;
  end;
  FuncCreateAccount:=function(EXT, GroupExt, ListObject)
    local FileNameFAC, FileNameEXT, FileNameINV, CompleteInfo, eInvariant, eInfoEXT, nbVert;
    eInfoEXT:=GroupFormalism.BankKeyInformation(EXT, GroupExt);
    nbVert:=Length(eInfoEXT.EXT);
    if nbRec=0 then
      MinNbVert:=nbVert;
    else
      MinNbVert:=Minimum(MinNbVert, nbVert);
    fi;
#    Print("We have eInfoEXT\n");
    eInvariant:=FuncInvariant(EXT);
#    Print("We have eInvariant\n");
    nbRec:=nbRec+1;
    #
    CompleteInfo:=GroupFormalism.BankCompleteInformation(EXT, GroupExt, ListObject);
    if DataBank.Saving then
      FileNameFAC:=Concatenation(DataBank.BankPath, "AccountFAC_", String(nbRec));
      SaveDataToFilePlusTouch(FileNameFAC, CompleteInfo);
    else
      Add(ListCompleteInformation, CompleteInfo);
    fi;
#    Print("After write FAC\n");
    #
    if DataBank.Saving then
      FileNameEXT:=Concatenation(DataBank.BankPath, "AccountEXT_", String(nbRec));
      SaveDataToFilePlusTouch(FileNameEXT, eInfoEXT);
    else
      Add(WRL, eInfoEXT);
    fi;
#    Print("After write EXT\n");
    #
    if DataBank.Saving then
      FileNameINV:=Concatenation(DataBank.BankPath, "AccountINV_", String(nbRec));
      SaveDataToFilePlusTouch(FileNameINV, eInvariant);
    fi;
    Add(ListInvariant, eInvariant);
#    Print("After write INV\n");
  end;
  FuncClearAccount:=function()
    local iRec, FileNameFAC, FileNameEXT, FileNameINV;
    if DataBank.Saving then
      for iRec in [1..nbRec]
      do
        FileNameFAC:=Concatenation(DataBank.BankPath, "AccountFAC_", String(iRec));
        FileNameEXT:=Concatenation(DataBank.BankPath, "AccountEXT_", String(iRec));
        FileNameINV:=Concatenation(DataBank.BankPath, "AccountINV_", String(iRec));
        RemoveFileIfExistPlusTouch(FileNameFAC);
        RemoveFileIfExistPlusTouch(FileNameEXT);
        RemoveFileIfExistPlusTouch(FileNameINV);
      od;
    else
      Unbind(WRL);
      Unbind(ListCompleteInformation);
    fi;
    Unbind(ListInvariant);
  end;
  return rec(FuncStabilizer:=FuncStabilizer, FuncCreateAccount:=FuncCreateAccount, FuncRetrieveObject:=FuncRetrieveObject, FuncClearAccount:=FuncClearAccount);
end;


#
# Data=rec(TheDepth:=...
#, ThePath:=..., IsBankSave:=false
#, IsBankSave:=...
#, IsRespawn:=...
#, GroupFormalism:=...
#, DualDescription:=...
#, Saving:=..., ThePathSave:=...
# This function computes the polyhedral description (given by incidences)
# of EXT and returns the orbits up to the symmetry group GivenSymmetry.
# GivenSymmetry should be used according to the group formalism specified in
#
# all tricks in the book are used:
# ---Bank formalism for storing used data
# ---Balinski theorem
# ---Recursive Adjacency decomposition method
# ---special symmetry groups of faces of the polytope
#
# saving system concerns only and exclusively the adjacency decomposition
# method itself, otherwise no save.
# the end result is that one can recover an interrupted computation
# with no error.
__ListFacetByAdjacencyDecompositionMethod:=function(EXT, GivenSymmetry, Data, BankFormalism)
  local RECListOrbit, IsFinished, eOrb, eInc, Ladj, SelectedOrbit, jOrb, MiniINCD, RPL, RPLift, testBank, OrdGRP, TheDim, WorkingSymGroup, LRES, NewData, RedStab, TheDate1, TheDate2, EllapsedTime, NewPathSave, TestNeedMoreSymmetry, ReturnedList, eSetUndone, nbUndone, testSym, fInc;
  if IsCorrectPath(Data.ThePathSave)=false then
    Print("Variable Data.ThePathSave=", Data.ThePathSave, " is incorrect\n");
    Error("It should finish with a /");
  fi;
  if IsDirectoryPath(Data.ThePathSave)=false and Data.Saving=true then
    Error("Directory Data.ThePathSave=", Data.ThePathSave, " is nonexistent");
  fi;
  if IsCorrectPath(Data.ThePath)=false then
    Print("Variable Data.ThePath=", Data.ThePath, " is incorrect\n");
    Error("It should finish with a /");
  fi;
  if IsDirectoryPath(Data.ThePath)=false then
    Error("Directory Data.ThePath=", Data.ThePath, " is nonexistent");
  fi;
#  Print("Before testBank\n");
  testBank:=BankFormalism.FuncRetrieveObject(EXT, GivenSymmetry);
#  Print("After testBank\n");
  if testBank<>false then
#    Print("Retrieve data from the bank\n");
    return Data.GroupFormalism.LiftingOrbits(EXT, testBank.ListOrbitFacet, GivenSymmetry, testBank.GRP);
  fi;
  TheDate1:=GetDate();
  # we would like to use IsBankSave but it is not possible with EllaspedTime
#  Print("Before TestNeedMoreSymmetry\n");
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
#  Print("After TestNeedMoreSymmetry\n");
  if testSym=true then
    WorkingSymGroup:=Data.GroupFormalism.GroupUnion(BankFormalism.FuncStabilizer(EXT), GivenSymmetry);
  else
    WorkingSymGroup:=GivenSymmetry;
  fi;
  OrdGRP:=Data.GroupFormalism.Order(WorkingSymGroup);
#  Print("OrdGRP=", OrdGRP, "\n");
  if Data.IsRespawn(OrdGRP, EXT, Data.TheDepth)=false then
    ReturnedList:=Data.DualDescriptionFunction(EXT, Data.GroupFormalism.ToPermGroup(EXT, GivenSymmetry), Data.ThePath);
    TheDate2:=GetDate();
    EllapsedTime:=TheDate2-TheDate1;
    if EllapsedTime > 10 then
      Print("EllapsedTime=", EllapsedTime, "\n");
    fi;
    if Data.IsBankSave(EllapsedTime, OrdGRP, EXT, Data.TheDepth)=true then
      Print("Before FuncCreateAccount\n");
      BankFormalism.FuncCreateAccount(EXT, GivenSymmetry, ReturnedList);
      Print("After FuncCreateAccount\n");
    fi;
  else
    TheDim:=RankMat(EXT)-1;
    Print("RESPAWN a new ADM computation |GRP|=", OrdGRP, " TheDim=", TheDim, " |EXT|=", Length(EXT), "\n");
#    FileSaveEXT:=Concatenation("EXT", String(Length(EXT)));
#    SaveDataToFile(FileSaveEXT, EXT);
    RPL:=Data.GroupFormalism.OrbitGroupFormalism(EXT, WorkingSymGroup, Data.ThePathSave, Data.Saving);
    Print("nbOrbit=", RPL.FuncNumberOrbit(), "\n");
    if RPL.FuncNumberOrbit()=0 then
      for eInc in Data.GetInitialRays(EXT, 10)
      do
        RPL.FuncInsert(eInc);
      od;
    fi;
    Print("Running the adjacency method recursively\n");
    while(true)
    do
      eSetUndone:=RPL.ComputeIntersectionUndone();
      nbUndone:=RPL.FuncNumberUndone();
      if RPL.FuncNumberOrbitDone()>0 then
        if nbUndone<=TheDim-1 or Length(eSetUndone)>0 then
          Print("End of computation, nbObj=", RPL.FuncNumber(), " nbOrbit=", RPL.FuncNumberOrbit(), " nbUndone=", nbUndone, " |eSetUndone|=", Length(eSetUndone), " Depth=", Data.TheDepth, " |EXT|=", Length(EXT), "\n");
          break;
        fi;
      fi;
      SelectedOrbit:=RPL.FuncGetMinimalUndoneOrbit();
      eInc:=RPL.FuncRecord(SelectedOrbit);
      Print("\n");
      RedStab:=Data.GroupFormalism.Stabilizer(EXT, WorkingSymGroup, eInc);
      Print("Considering orbit ", SelectedOrbit, " |inc|=", Length(eInc), " Depth=", Data.TheDepth, " |stab|=", Order(RedStab), " dim=", TheDim, "\n");
      RPLift:=__ProjectionLiftingFramework(EXT, eInc);
      NewPathSave:=Concatenation(Data.ThePathSave, "OrbitRespawn", String(SelectedOrbit), "/");
      CreateDirectoryPlusTest(NewPathSave, Data.Saving);
      NewData:=ShallowCopy(Data);
      NewData.TheDepth:=NewData.TheDepth+1;
      NewData.ThePathSave:=NewPathSave;
      Ladj:=__ListFacetByAdjacencyDecompositionMethod(EXT{eInc}, RedStab, NewData, BankFormalism);
      Print("We treat ", Length(Ladj), " orbits\n");
      RemoveDirectoryPlusTest(NewPathSave, Data.Saving);
      for fInc in Ladj
      do
        RPL.FuncInsert(RPLift.FuncLift(fInc));
      od;
      RPL.FuncPutOrbitAsDone(SelectedOrbit);
      nbUndone:=RPL.FuncNumberUndone();
      Print("We have ", RPL.FuncNumberOrbit(), " orbits");
      Print("  Nr treated=", RPL.FuncNumberOrbitDone(), " orbits");
      Print("  nbUndone=", nbUndone);
      Print("\n");
    od;
    LRES:=RPL.FuncListOrbitIncidence();
    ReturnedList:=Data.GroupFormalism.LiftingOrbits(EXT, LRES, GivenSymmetry, WorkingSymGroup);
    Print("LRES=", Length(LRES), " |ReturnedList|=", Length(ReturnedList), " |GivenSymmetry|=", Order(GivenSymmetry), " |WorkingSymGroup|=", Order(WorkingSymGroup), "\n");
    TheDate2:=GetDate();
    EllapsedTime:=TheDate2-TheDate1;
    if EllapsedTime > 10 then
      Print("EllapsedTime=", EllapsedTime, "\n");
    fi;
    if Data.IsBankSave(EllapsedTime, OrdGRP, EXT, Data.TheDepth)=true then
      Print("Before FuncCreateAccount\n");
      BankFormalism.FuncCreateAccount(EXT, WorkingSymGroup, LRES);
      Print("After FuncCreateAccount\n");
    fi;
  fi;
  return ReturnedList;
end;
