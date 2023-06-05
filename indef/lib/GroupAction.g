

ReducePrimeMultiplicity:=function(TheValue)
    local ListInt, SetInt;
    if TheValue = 1 then
        return 1;
    fi;
    ListInt:=FactorsInt(TheValue);
    SetInt:=Set(ListInt);
    return Product(SetInt);
end;


GetRationalInvariant:=function(GRPmatr)
    local ListPrimes, eGen, TheDen, ListInt, SetInt;
    ListPrimes:=[];
    for eGen in GeneratorsOfGroup(GRPmatr)
    do
        if IsIntegralMat(eGen)=false then
            TheDen:=GetDenominatorMatrix(eGen);
            ListInt:=FactorsInt(TheDen);
            SetInt:=Set(ListInt);
            Append(ListPrimes, SetInt);
        fi;
    od;
    if Length(ListPrimes) = 0 then
        return 1;
    fi;
    return Product(Set(ListPrimes));
end;


LLLMatrixGroupReduction:=function(n, GRPmatr)
    local PosDefMat, LGen, eMat, eProd, eRec, Pmat, PmatInv, ListMatrNew, eMatNew;
    if IndefinitePrint then
        Print("Begin of LLLMatrixGroupReduction\n");
    fi;
    PosDefMat:=IdentityMat(n);
    LGen:=GeneratorsOfGroup(GRPmatr);
    for eMat in LGen
    do
        if eMat<>IdentityMat(n) then
            eProd:=eMat * TransposedMat(eMat);
            PosDefMat:=PosDefMat + eProd;
        fi;
    od;
    eRec:=LLLReducedGramMat(PosDefMat);
    Pmat:=eRec.transformation;
    PmatInv:=Inverse(Pmat);
    ListMatrNew:=[];
    for eMat in LGen
    do
        eMatNew:=Pmat * eMat * PmatInv;
        Add(ListMatrNew, eMatNew);
    od;
    if IndefinitePrint then
        Print("End of LLLMatrixGroupReduction\n");
    fi;
    return rec(GRPred:=Group(ListMatrNew), Pmat:=Pmat);
end;


CosetRepresentative_Stabilizer_TwoAct:=function(Id1, ListGen1, ListGen2, ePt, TheAct)
    local ListPt, ListCoset, nbGen, FuncAction, ThePos, LastPos, pos, aPt, iGen, eGen1, eGen2, uPt;
    ListPt:=[ePt];
    ListCoset:=[Id1];
    nbGen:=Length(ListGen1);
    FuncAction:=function(uPt, eGen1)
        if Position(ListPt, uPt)=fail then
            Add(ListPt, uPt);
            Add(ListCoset, eGen1);
        fi;
    end;
    ThePos:=0;
    LastPos:=1;
    while(true)
    do
        for pos in [ThePos+1..LastPos]
        do
            aPt:=ListPt[pos];
            for iGen in [1..nbGen]
            do
                eGen1:=ListCoset[pos] * ListGen1[iGen];
                eGen2:=ListGen2[iGen];
                uPt:=TheAct(aPt, eGen2);
                FuncAction(uPt, eGen1);
            od;
        od;
        ThePos:=LastPos;
        LastPos:=Length(ListPt);
        if ThePos=LastPos then
            break;
        fi;
    od;
    return ListCoset;
end;


LinearSpace_GetDivisor:=function(TheSpace)
  local n, TheDet, eDiv, eMat, IsOK, eVect, eSol;
  n:=Length(TheSpace);
  TheDet:=AbsInt(DeterminantMat(TheSpace));
  eDiv:=1;
  while(true)
  do
    eMat:=eDiv*IdentityMat(n);
    IsOK:=true;
    for eVect in eMat
    do
      eSol:=SolutionIntMat(TheSpace, eVect);
      if eSol=fail then
        IsOK:=false;
      fi;
    od;
    if IsOK=true then
      return eDiv;
    fi;
    if IndefinitePrint then
        Print("eDiv=", eDiv, "\n");
    fi;
    eDiv:=eDiv+1;
  od;
  Error("We should never reach that stage");
end;


OrbitComputationGeneral_limited:=function(GRPmatr, start_obj, TheAction, MaxSize)
    local TheList, TheDict, f_insert, ListGen, CurrPos, len, eGen, idx, new_obj;
    TheList:=[];
    TheDict:=NewDictionary(start_obj, true);
    f_insert:=function(eVect)
        if LookupDictionary(TheDict, eVect)<>fail then
            return;
        fi;
        Add(TheList, eVect);
        AddDictionary(TheDict, eVect, true);
    end;
    f_insert(start_obj);
    ListGen:=GeneratorsOfGroup(GRPmatr);
    CurrPos:=1;
    while(true)
    do
        len:=Length(TheList);
        if IndefinitePrint then
            Print("OrbitComputationGeneral_limited len=", len, "\n");
        fi;
        for eGen in ListGen
        do
            for idx in [CurrPos..len]
            do
                new_obj:=TheAction(TheList[idx], eGen);
                f_insert(new_obj);
            od;
            if MaxSize > 0 and Length(TheList) > MaxSize then
                if IndefinitePrint then
                    Print("Orbit enumeration terminated because of too large size\n");
                fi;
                return fail;
            fi;
        od;
        if Length(TheList)=len then
            break;
        fi;
        CurrPos:=len+1;
    od;
    return Set(TheList);
end;


OrbitComputation_limited:=function(GRPmatr, start_vect, TheMod, MaxSize)
    local TheList, TheDict, f_insert, ListGen, CurrPos, len, eGen, idx, eVect1, eVect2;
    TheList:=[];
    TheDict:=NewDictionary(start_vect, true);
    f_insert:=function(eVect)
        if LookupDictionary(TheDict, eVect)<>fail then
            return;
        fi;
        Add(TheList, eVect);
        AddDictionary(TheDict, eVect, true);
    end;
    f_insert(start_vect);
    ListGen:=GeneratorsOfGroup(GRPmatr);
    CurrPos:=1;
    while(true)
    do
        len:=Length(TheList);
        if IndefinitePrint then
            Print("OrbitComputation_limited len=", len, "\n");
        fi;
        for eGen in ListGen
        do
            for idx in [CurrPos..len]
            do
                eVect1:=TheList[idx] * eGen;
                eVect2:=List(eVect1, x->x mod TheMod);
                f_insert(eVect2);
            od;
            if MaxSize > 0 and Length(TheList) > MaxSize then
                if IndefinitePrint then
                    Print("Orbit enumeration terminated because of too large size\n");
                fi;
                return fail;
            fi;
        od;
        if Length(TheList)=len then
            break;
        fi;
        CurrPos:=len+1;
    od;
    return Set(TheList);
end;


IsStabilizingMod:=function(TheGRP, RecSpace)
  local eVect, eGen, eSol;
  for eVect in RecSpace.TheSpace
  do
    for eGen in GeneratorsOfGroup(TheGRP)
    do
      eSol:=SolutionIntMat(RecSpace.TheSpaceMod, eVect*eGen);
      if eSol=fail then
        return List(eVect, x->x mod RecSpace.TheMod);
      fi;
    od;
  od;
  return true;
end;


IsStabilizing:=function(TheGRP, TheSpace)
  local eVect, eGen, eSol;
  for eVect in TheSpace
  do
    for eGen in GeneratorsOfGroup(TheGRP)
    do
      eSol:=SolutionIntMat(TheSpace, eVect*eGen);
      if eSol=fail then
        return false;
      fi;
    od;
  od;
  return true;
end;


#GetListPermGenOrderedOrbit
MapToPermutationOrderedOrbit:=function(eInput, TheAction, O)
  local GetListPermGenOrderedOrbit;
  GetListPermGenOrderedOrbit:=function(ListGens)
    local ListPermGens, eGen, eListImg, ePerm1, eList, ePerm2;
    ListPermGens:=[];
    for eGen in ListGens
    do
      eListImg:=List(O, x->TheAction(x, eGen));
      ePerm1:=SortingPerm(eListImg);
#      eList:=List(O, x->Position(O, TheAction(x, eGen)));
#      ePerm2:=PermList(eList);
#      if ePerm1<>ePerm2 then
#        Error("ePerm1 <> ePerm2");
#      fi;
      Add(ListPermGens, ePerm1);
    od;
    return ListPermGens;
  end;
  if IsList(eInput) then
    return GetListPermGenOrderedOrbit(eInput);
  fi;
  if IsGroup(eInput) then
    return PersoGroupPerm(GetListPermGenOrderedOrbit(GeneratorsOfGroup(eInput)));
  fi;
  Error("Failed to find a mathing type");
end;


LinearSpace_ModStabilizer:=function(GRPmatr, TheSpace, TheMod)
    local n, TheSpaceMod, RecSpace, TheAction, GRPret, test, O, ListMatrGens, ListPermGens, eGen, eList, GRPperm, eSet, eStab, phi, ePerm1, ePerm2, eListImg, MaxSize;
    n:=Length(TheSpace);
    if IndefinitePrint then
        Print("TheMod=", TheMod, "\n");
    fi;
    TheSpaceMod:=Concatenation(TheSpace, TheMod*IdentityMat(n));
    RecSpace:=rec(TheSpace:=TheSpace, TheSpaceMod:=TheSpaceMod, TheMod:=TheMod);
    TheAction:=function(eClass, eElt)
        local eVect;
        eVect:=eClass*eElt;
        return List(eVect, x->x mod TheMod);
    end;
    GRPret:=Group(GeneratorsOfGroup(GRPmatr));
    MaxSize:=100000;
    if IndefinitePrint then
        Print("Order(GRPmatr)=", Order(GRPmatr), "\n");
    fi;
    while(true)
    do
        test:=IsStabilizingMod(GRPret, RecSpace);
        if test=true then
            return GRPret;
        fi;
        O:=OrbitComputation_limited(GRPret, test, TheMod, MaxSize);
        if O=fail then
            return -1;
        fi;
        ListMatrGens:=GeneratorsOfGroup(GRPret);
        ListPermGens:=MapToPermutationOrderedOrbit(ListMatrGens, TheAction, O);
        eSet:=Filtered([1..Length(O)], x->SolutionIntMat(TheSpaceMod, O[x])<>fail);
        if IndefinitePrint then
            Print("|O|=", Length(O), " |eSet|=", Length(eSet), "\n");
        fi;
        GRPret:=Stabilizer(GRPret, eSet, ListMatrGens, ListPermGens, OnSets);
        if IndefinitePrint then
            Print("|O|=", Length(O), " |eSet|=", Length(eSet), " |ListPermGens|=", Length(ListPermGens), " |ListGen(GRPret)|=", Length(GeneratorsOfGroup(GRPret)), " |GRPperm|=", Order(Group(ListPermGens)), " |Oset|=", Length(Orbit(Group(ListPermGens), eSet, OnSets)), "\n");
        fi;
    od;
end;


LinearSpace_ModStabilizer_RightCoset:=function(RecMatr, TheSpace, TheMod)
  local n, TheSpaceMod, RecSpace, TheAction, GRPret, test, O, ListMatrGens, ListPermGens, eGen, eList, GRPperm, eSet, eStab, ListListCoset, Id1, ListCoset, ePerm1, ePerm2, eListImg;
  n:=Length(TheSpace);
  TheSpaceMod:=Concatenation(TheSpace, TheMod*IdentityMat(n));
  RecSpace:=rec(TheSpace:=TheSpace, TheSpaceMod:=TheSpaceMod, TheMod:=TheMod);
  TheAction:=function(eClass, eElt)
    local eVect;
    eVect:=eClass*eElt;
    return List(eVect, x->x mod TheMod);
  end;
  GRPret:=PersoGroup(GeneratorsOfGroup(RecMatr.GRPmatr), Identity(RecMatr.GRPmatr));
  ListListCoset:=ShallowCopy(RecMatr.ListListCoset);
  Id1:=IdentityMat(n);
  while(true)
  do
    test:=IsStabilizingMod(GRPret, RecSpace);
    if test=true then
      return rec(ListListCoset:=ListListCoset, GRPmatr:=GRPret);
    fi;
    O:=Set(Orbit(GRPret, test, TheAction));
    if IndefinitePrint then
        Print("|O|=", Length(O), "\n");
    fi;
    ListMatrGens:=GeneratorsOfGroup(GRPret);
    ListPermGens:=MapToPermutationOrderedOrbit(ListMatrGens, TheAction, O);
    eSet:=Filtered([1..Length(O)], x->SolutionIntMat(TheSpaceMod, O[x])<>fail);
    GRPret:=Stabilizer(GRPret, eSet, ListMatrGens, ListPermGens, OnSets);
    ListCoset:=CosetRepresentative_Stabilizer_TwoAct(Id1, ListMatrGens, ListPermGens, eSet, OnSets);
    Add(ListListCoset, ListCoset);
  od;
end;


LinearSpace_Stabilizer_RightCoset:=function(GRPmatr, TheSpace_pre)
    local TheSpace, LFact, eList, GRPret, TheMod, i, RecGRP;
    TheSpace:=LLLReducedBasis(TheSpace_pre).basis;
    RecGRP:=rec(ListListCoset:=[], GRPmatr:=PersoGroup(GeneratorsOfGroup(GRPmatr), Identity(GRPmatr)));
    if IsStabilizing(GRPmatr, TheSpace) then
        return RecGRP;
    fi;
    LFact:=LinearSpace_GetDivisor(TheSpace);
    eList:=FactorsInt(LFact);
    if IndefinitePrint then
        Print("LFact=", LFact, " eList=", eList, "\n");
    fi;
    for i in [1..Length(eList)]
    do
        TheMod:=Product(eList{[1..i]});
        RecGRP:=LinearSpace_ModStabilizer_RightCoset(RecGRP, TheSpace, TheMod);
        if IsStabilizing(RecGRP.GRPmatr, TheSpace) then
            return RecGRP;
        fi;
    od;
    if IsStabilizing(RecGRP.GRPmatr, TheSpace)=false then
        Error("Algorithm error");
    fi;
    return RecGRP;
end;


LinearSpace_ExpandListListCoset:=function(n, ListListCoset)
    local ListCosetRet, eListCoset, ListCosetNew, eCos, fCos;
    ListCosetRet:=[IdentityMat(n)];
    for eListCoset in ListListCoset
    do
        ListCosetNew:=[];
        for eCos in ListCosetRet
        do
            for fCos in eListCoset
            do
                Add(ListCosetNew, eCos * fCos);
            od;
        od;
        ListCosetRet:=ListCosetNew;
    od;
    return ListCosetRet;
end;


LinearSpace_Stabilizer_Direct:=function(GRPmatr, TheSpace)
    local TheSpaceCan, TheAction;
    if IndefinitePrint then
        Print("Beginning of LinearSpace_Stabilizer_Direct\n");
    fi;
    TheSpaceCan:=HermiteNormalFormIntegerMat(TheSpace);
    TheAction:=function(eSpace, eElt)
        return HermiteNormalFormIntegerMat(eSpace * eElt);
    end;
    if IndefinitePrint then
        Print("|O|=", Length(OrbitComputationGeneral_limited(GRPmatr, TheSpaceCan, TheAction, -1)), "\n");
    fi;
    return Stabilizer(GRPmatr, TheSpaceCan, TheAction);
end;


LinearSpace_Equivalence_Direct:=function(GRPmatr, eSpace, fSpace)
    local eSpaceCan, fSpaceCan, TheAction, TheEquiv;
    eSpaceCan:=HermiteNormalFormIntegerMat(eSpace);
    fSpaceCan:=HermiteNormalFormIntegerMat(fSpace);
    TheAction:=function(uSpace, eElt)
        return HermiteNormalFormIntegerMat(uSpace * eElt);
    end;
    TheEquiv:=RepresentativeAction(GRPmatr, eSpaceCan, fSpaceCan, TheAction);
    return TheEquiv;
end;


#
#
# for L a linear space of finite index in Z^n
LinearSpace_Stabilizer_Kernel_Reduced:=function(GRPmatr, TheSpace_pre)
  local TheSpace, LFact, eList, GRPret, TheMod, i, result;
  TheSpace:=LLLReducedBasis(TheSpace_pre).basis;
  if IsStabilizing(GRPmatr, TheSpace) then
    return GRPmatr;
  fi;
  GRPret:=PersoGroup(GeneratorsOfGroup(GRPmatr), Identity(GRPmatr));
  LFact:=LinearSpace_GetDivisor(TheSpace);
  if IndefinitePrint then
      Print("LFact=", LFact, "\n");
  fi;
  eList:=FactorsInt(LFact);
  for i in [1..Length(eList)]
  do
    TheMod:=Product(eList{[1..i]});
    result:=LinearSpace_ModStabilizer(GRPret, TheSpace, TheMod);
    if result=-1 then
        return LinearSpace_Stabilizer_Direct(GRPret, TheSpace);
    fi;
    GRPret := result;
    if IsStabilizing(GRPret, TheSpace) then
      return GRPret;
    fi;
  od;
  if IsStabilizing(GRPret, TheSpace)=false then
    Error("Algorithm error");
  fi;
  return GRPret;
end;


LinearSpace_ComputeOrbit:=function(GRPmatr, TheSpace)
    local ListSpaces, n_new, n_old, f_equal, f_insert, ListGen, CurrentPos, len, i, eGen, NewSpace, NewLen;
    ListSpaces:=[];
    n_new:=0;
    n_old:=0;
    f_insert:=function(eSpace)
        local eSpaceCan, quot;
        # SpaceCan = P * eSpace
        eSpaceCan:=HermiteNormalFormIntegerMat(eSpace);
        if Position(ListSpaces, eSpaceCan)<>fail then
            n_old:=n_old+1;
            return;
        fi;
        n_new:=n_new+1;
        Add(ListSpaces, eSpaceCan);
        quot:=(n_old + 0.0) / (n_new + 0.0);
        if IndefinitePrint then
            Print("Now |ListSpaces|=", Length(ListSpaces), " n_old=", n_old, " n_new=", n_new, " quot=", quot, "\n");
        fi;
    end;
    f_insert(TheSpace);
    ListGen:=GeneratorsOfGroup(GRPmatr);
    CurrentPos:=1;
    while(true)
    do
        len:=Length(ListSpaces);
        for i in [CurrentPos..len]
        do
            for eGen in ListGen
            do
                NewSpace:=ListSpaces[i] * eGen;
                f_insert(NewSpace);
            od;
        od;
        NewLen:=Length(ListSpaces);
        if len = NewLen then
            break;
        fi;
        CurrentPos:=len+1;
    od;
    return ListSpaces;
end;


LinearSpace_Stabilizer_Kernel:=function(GRPmatr, TheSpace_pre)
    local n, eRec, Pmat, PmatInv, TheSpace_B, TheSpace_C, GRPret, ListGenNew, eGen, eGenNew;
    n:=Length(TheSpace_pre);
    eRec:=LLLMatrixGroupReduction(n, GRPmatr);
    Pmat:=eRec.Pmat;
    PmatInv:=Inverse(Pmat);
    TheSpace_B:=TheSpace_pre * PmatInv;
    TheSpace_C:=LLLbasisReduction(TheSpace_B).LattRed;
    if IndefinitePrint then
        Print("|ListSpaces|=", Length(LinearSpace_ComputeOrbit(eRec.GRPred, TheSpace_C)), "\n");
    fi;
    GRPret:=LinearSpace_Stabilizer_Kernel_Reduced(eRec.GRPred, TheSpace_C);
    ListGenNew:=[];
    for eGen in GeneratorsOfGroup(GRPret)
    do
        eGenNew:=PmatInv * eGen * Pmat;
        Add(ListGenNew, eGenNew);
    od;
    return Group(ListGenNew);
end;


LinearSpace_Stabilizer:=function(GRPmatr, TheSpace_pre)
    local result, TheAnswer;
    result:=LinearSpace_Stabilizer_Kernel(GRPmatr, TheSpace_pre);
    TheAnswer:=rec(GRPmatr:=GRPmatr, TheSpace_pre:=TheSpace_pre, result:=result);
    return result;
end;


LinearSpace_ModEquivalence:=function(GRPmatr, NeedStabilizer, TheSpace1, TheSpace2, TheMod)
    local MaxSize, n, TheSpace1Mod, TheSpace2Mod, RecSpace2, TheAction, IsEquiv, GRPwork, eElt, test, eVect, O, ListMatrGens, ListPermGens, eGen, eList, GRPperm, eSet1, eSet2, eTest, eStab, eMat, TheSpace1work, GenerateGroupInfo, GetFace, GrpInf, test1, test2;
    if IndefinitePrint then
        Print("LinearSpace_ModEquivalence, TheMod=", TheMod, "\n");
        Print("TheSpace1=\n");
        PrintArray(TheSpace1);
        Print("TheSpace2=\n");
        PrintArray(TheSpace2);
        Print("det(TheSpace1)=", DeterminantMat(TheSpace1), " det(TheSpace2)=", DeterminantMat(TheSpace2), "\n");
    fi;
    MaxSize:=100000;
    n:=Length(TheSpace1);
    TheSpace1Mod:=Concatenation(TheSpace1, TheMod*IdentityMat(n));
    TheSpace2Mod:=Concatenation(TheSpace2, TheMod*IdentityMat(n));
    RecSpace2:=rec(TheSpace:=TheSpace2, TheSpaceMod:=TheSpace2Mod, TheMod:=TheMod);
    TheAction:=function(eClass, eElt)
        local eVect;
        eVect:=eClass*eElt;
        return List(eVect, x->x mod TheMod);
    end;
    IsEquiv:=function(eEquiv)
        local eVect, eGen, eSol;
        for eVect in TheSpace1
        do
            eSol:=SolutionIntMat(TheSpace2Mod, eVect*eEquiv);
            if eSol=fail then
                return List(eVect*eEquiv, x->x mod TheMod);
            fi;
        od;
        return true;
    end;
    GenerateGroupInfo:=function(test)
        local ListMatrGens, ListPermGens, eGen, eList, eListImg, ePerm1, ePerm2;
        O:=OrbitComputation_limited(GRPwork, test, TheMod, MaxSize);
        if O=fail then
            return -1;
        fi;
        if IndefinitePrint then
            Print("|O|=", Length(O), "\n");
        fi;
        ListMatrGens:=GeneratorsOfGroup(GRPwork);
        ListPermGens:=[];
        for eGen in ListMatrGens
        do
            eListImg:=List(O, x->TheAction(x, eGen));
            ePerm1:=SortingPerm(eListImg);
            eList:=List(O, x->Position(O, TheAction(x, eGen)));
            ePerm2:=PermList(eList);
            if ePerm1<>ePerm2 then
                Error("ePerm1 <> ePerm2");
            fi;
            Add(ListPermGens, ePerm2);
        od;
        if IndefinitePrint then
            Print("ListPermGens built\n");
        fi;
        return rec(O:=O, ListMatrGens:=ListMatrGens, ListPermGens:=ListPermGens);
    end;
    GetFace:=function(TheO, TheSpace)
        return Filtered([1..Length(TheO)], x->SolutionIntMat(TheSpace, TheO[x])<>fail);
    end;
    GRPwork:=PersoGroup(GeneratorsOfGroup(GRPmatr), Identity(GRPmatr));
    eElt:=IdentityMat(n);
    while(true)
    do
        if IndefinitePrint then
            Print("Before test1, test2 equivalence\n");
        fi;
        test1:=IsEquiv(eElt);
        if NeedStabilizer=false then
            test2:=true;
        else
            test2:=IsStabilizingMod(GRPwork, RecSpace2);
        fi;
        if IndefinitePrint then
            Print("test1=true=", test1=true, " test2=true=", test2=true, "\n");
        fi;
        if test1=true and test2=true then
            if IndefinitePrint then
                Print("Returning from LinearSpace_ModEquivalence\n");
            fi;
            return rec(GRPwork:=GRPwork, eEquiv:=eElt);
        fi;
        if test1<>true then
            GrpInf:=GenerateGroupInfo(test1);
            if GrpInf=-1 then
                return -1;
            fi;
            TheSpace1work:=TheSpace1Mod*eElt;
            eSet1:=GetFace(GrpInf.O, TheSpace1work);
            eSet2:=GetFace(GrpInf.O, TheSpace2Mod);
            if IndefinitePrint then
                Print("|eSet1|=", Length(eSet1), " |eSet2|=", Length(eSet2), " |GrpInf.O|=", Length(GrpInf.O), " |GrpInf.ListMatrGens|=", Length(GrpInf.ListMatrGens), "\n");
            fi;
            eMat:=RepresentativeAction(GRPwork, eSet1, eSet2, GrpInf.ListMatrGens, GrpInf.ListPermGens, OnSets);
            if IndefinitePrint then
                Print("We have eMat\n");
            fi;
            if eMat=fail then
                return fail;
            fi;
            if IndefinitePrint then
                Print("Before computing GRPwork\n");
            fi;
            GRPwork:=Stabilizer(GRPwork, eSet2, GrpInf.ListMatrGens, GrpInf.ListPermGens, OnSets);
            if IndefinitePrint then
                Print("After stabilization |GRPwork|=", Order(GRPwork), "\n");
            fi;
            eElt:=eElt*eMat;
        fi;
        if test2<>true then
            if IndefinitePrint then
                Print("Before GrpInf\n");
            fi;
            GrpInf:=GenerateGroupInfo(test2);
            if GrpInf=-1 then
                return -1;
            fi;
            if IndefinitePrint then
                Print("We have GrpInf\n");
            fi;
            eSet2:=GetFace(GrpInf.O, TheSpace2Mod);
            if IndefinitePrint then
                Print("|eSet2|=", Length(eSet2), " |GrpInf.O|=", Length(GrpInf.O), "\n");
            fi;
            GRPwork:=Stabilizer(GRPwork, eSet2, GrpInf.ListMatrGens, GrpInf.ListPermGens, OnSets);
            if IndefinitePrint then
                Print("We have GRPwork\n");
            fi;
        fi;
    od;
end;


LinearSpace_Equivalence_Kernel_Reduced:=function(GRPmatr, TheSpace1_pre, TheSpace2_pre)
    local TheSpace1, TheSpace2, n, LFact1, LFact2, eList, IsEquivalence, GRPwork, eElt, TheMod, TheSpace1Img, eTest, i, eDet1, eDet2, NeedStabilizer;
    TheSpace1:=LLLReducedBasis(TheSpace1_pre).basis;
    TheSpace2:=LLLReducedBasis(TheSpace2_pre).basis;
    n:=Length(TheSpace1);
    LFact1:=LinearSpace_GetDivisor(TheSpace1);
    LFact2:=LinearSpace_GetDivisor(TheSpace2);
    if IndefinitePrint then
        eDet1:=AbsInt(DeterminantMat(TheSpace1));
        eDet2:=AbsInt(DeterminantMat(TheSpace2));
        Print("eDet1=", eDet1, " eDet2=", eDet2, "\n");
        Print("LFact1=", LFact1, " LFact2=", LFact2, "\n");
    fi;
    if LFact1<>LFact2 then
        return fail;
    fi;
    eList:=FactorsInt(LFact1);
    IsEquivalence:=function(eEquiv)
        local eVect, eSol;
        if IndefinitePrint then
            Print("Beginning of IsEquivalence\n");
        fi;
        for eVect in TheSpace1
        do
            eSol:=SolutionIntMat(TheSpace2, eVect*eEquiv);
            if eSol=fail then
                if IndefinitePrint then
                    Print("Returning false from IsEquivalence\n");
                fi;
                return false;
            fi;
        od;
        if IndefinitePrint then
            Print("Returning true from IsEquivalence\n");
        fi;
        return true;
    end;
    GRPwork:=PersoGroup(GeneratorsOfGroup(GRPmatr), Identity(GRPmatr));
    eElt:=IdentityMat(n);
    for i in [1..Length(eList)]
    do
        TheMod:=Product(eList{[1..i]});
        if IndefinitePrint then
            Print("i=", i, " TheMod=", TheMod, "\n");
        fi;
        if IsEquivalence(eElt) then
            if IndefinitePrint then
                Print("Returning eElt 1\n");
            fi;
            return eElt;
        fi;
        TheSpace1Img:=List(TheSpace1, x->x*eElt);
        NeedStabilizer:=true;
        if i = Length(eList) then
            NeedStabilizer:=false;
        fi;
        eTest:=LinearSpace_ModEquivalence(GRPwork, NeedStabilizer, TheSpace1Img, TheSpace2, TheMod);
        if eTest=-1 then
            if IndefinitePrint then
                Print("Using the direct approach\n");
            fi;
            eTest:=LinearSpace_Equivalence_Direct(GRPwork, TheSpace1Img, TheSpace2);
            if eTest=fail then
                if IndefinitePrint then
                    Print("Returning fail by direct approach\n");
                fi;
                return fail;
            fi;
            eElt:=eElt * eTest;
            if IsEquivalence(eElt)=false then
                Error("Algorithm error 1");
            fi;
            return eElt;
        fi;
        if IndefinitePrint then
            Print("We have eTest\n");
        fi;
        if eTest=fail then
            if IndefinitePrint then
                Print("Returning fail\n");
            fi;
            return fail;
        fi;
        eElt:=eElt*eTest.eEquiv;
        if NeedStabilizer then
            GRPwork:=eTest.GRPwork;
        fi;
    od;
    if IsEquivalence(eElt)=false then
        Error("Algorithm error 2");
    fi;
    if IndefinitePrint then
        Print("Returning eElt 2\n");
    fi;
    return eElt;
end;


LinearSpace_Equivalence_Kernel:=function(GRPmatr, TheSpace1_pre, TheSpace2_pre)
    local n, eRec, Pmat, PmatInv, TheSpace1_B, TheSpace2_B, TheSpace1_C, TheSpace2_C, opt, RetMat;
    n:=Length(TheSpace1_pre);
    if IndefinitePrint then
        Print("Generators(GRPmatr)=", GeneratorsOfGroup(GRPmatr), "\n");
    fi;
    eRec:=LLLMatrixGroupReduction(n, GRPmatr);
    if IndefinitePrint then
        Print("Generators(eRec.GRPred)=", GeneratorsOfGroup(eRec.GRPred), "\n");
    fi;
    Pmat:=eRec.Pmat;
    PmatInv:=Inverse(Pmat);
    TheSpace1_B:=TheSpace1_pre * PmatInv;
    TheSpace2_B:=TheSpace2_pre * PmatInv;
    TheSpace1_C:=LLLbasisReduction(TheSpace1_B).LattRed;
    TheSpace2_C:=LLLbasisReduction(TheSpace2_B).LattRed;
    if IndefinitePrint then
        Print("TheSpace1_C=", TheSpace1_C, "\n");
        Print("TheSpace2_C=", TheSpace2_C, "\n");
    fi;
    opt:=LinearSpace_Equivalence_Kernel_Reduced(eRec.GRPred, TheSpace1_C, TheSpace2_C);
    if opt=fail then
        return fail;
    fi;
    RetMat:=PmatInv * opt * Pmat;
    return RetMat;
end;


LinearSpace_Equivalence:=function(GRPmatr, TheSpace1_pre, TheSpace2_pre)
    local result, TheAnswer;
    result:=LinearSpace_Equivalence_Kernel(GRPmatr, TheSpace1_pre, TheSpace2_pre);
    TheAnswer:=rec(GRPmatr:=GRPmatr, TheSpace1_pre:=TheSpace1_pre, TheSpace2_pre:=TheSpace2_pre, result:=result);
    return result;
end;


# We want to find an invariant lattice for the group.
# This allows to conjugate the matrix group into an integral group.
#
# Such a lattice do not necessarily exists:
# ---For example if there are matrices of determinant a with |a| > 1.
# then no lattice exist.
#
# Such a lattice exist in many cases:
# ---For example if the group is finite.
# ---In application case such as indefinite forms and
#
#
MatrixIntegral_GetInvariantSpace:=function(n, GRPrat)
    local LGen, LGenTot, TheSpace, TheDet, IncreaseSpace;
    LGen:=GeneratorsOfGroup(GRPrat);
    LGenTot:=Set(Concatenation(LGen, List(LGen, Inverse)));
    if IndefinitePrint then
        Print("LGen|=", Length(LGen), " |LGenTot|=", Length(LGenTot), "\n");
    fi;
    TheSpace:=IdentityMat(n);
    TheDet:=1;
    IncreaseSpace:=function()
        local eGen, ConcatSpace, NewSpace, NewDet;
        for eGen in LGenTot
        do
            ConcatSpace:=Concatenation(TheSpace, TheSpace * eGen);
            NewSpace:=GetZbasis(ConcatSpace);
            NewDet:=AbsInt(DeterminantMat(NewSpace));
            if NewDet<>TheDet then
                TheSpace:=ShallowCopy(NewSpace);
                TheDet:=NewDet;
                return false;
            fi;
        od;
        return true;
    end;
    while(true)
    do
        if IncreaseSpace() then
            break;
        fi;
    od;
    return TheSpace;
end;


InvariantSpaceOfGroup:=function(n, GRPmatr)
    local LGen, LSeq, EquaSyst, NSP, eVect, eGen;
    LGen:=GeneratorsOfGroup(GRPmatr);
    if Length(LGen)=0 then
        return IdentityMat(n);
    fi;
    LSeq:=List(LGen, x->TransposedMat(x) - IdentityMat(n));
    EquaSyst:=TransposedMat(Concatenation(LSeq));
    NSP:=NullspaceMat(EquaSyst);
    for eVect in NSP
    do
        for eGen in LGen
        do
            if eVect * eGen<>eVect then
                Error("Vector is not invariant");
            fi;
        od;
    od;
    return NSP;
end;


MatrixIntegral_Stabilizer:=function(n, GRPrat)
    local LGen, TheSpace, TheSpaceInv, ListGenInt, GRPint, eStab, ListGenIntSpace;
    LGen:=GeneratorsOfGroup(GRPrat);
    TheSpace:=MatrixIntegral_GetInvariantSpace(n, GRPrat);
    TheSpaceInv:=Inverse(TheSpace); # Also known as LattToStab
    ListGenInt:=List(LGen, x->TheSpace * x * TheSpaceInv);
    GRPint:=Group(ListGenInt);
    if First(ListGenInt, x->IsIntegralMat(x)=false)<>fail then
        Error("Some geneatorsmatrix are not integral");
    fi;
    if IsIntegralMat(TheSpaceInv)=false then
        Error("TheSpaceInv is not integral");
    fi;
    eStab:=LinearSpace_Stabilizer(GRPint, TheSpaceInv);
    ListGenIntSpace:=List(GeneratorsOfGroup(eStab), x->TheSpaceInv * x * TheSpace);
    return PersoGroupMatrix(ListGenIntSpace, n);
end;


MatrixIntegral_RightCosets:=function(n, GRPrat)
    local LGen, TheSpace, TheSpaceInv, ListGenInt, GRPint, RecStab_RightCoset, ListCoset, ListCosetRet;
    LGen:=GeneratorsOfGroup(GRPrat);
    TheSpace:=MatrixIntegral_GetInvariantSpace(n, GRPrat);
    TheSpaceInv:=Inverse(TheSpace);
    ListGenInt:=List(LGen, x->TheSpace * x * TheSpaceInv);
    GRPint:=Group(ListGenInt);
    RecStab_RightCoset:=LinearSpace_Stabilizer_RightCoset(GRPint, TheSpaceInv);
    if IndefinitePrint then
        Print("RecStab_RightCoset=", RecStab_RightCoset, "\n");
    fi;
    ListCoset:=LinearSpace_ExpandListListCoset(n, RecStab_RightCoset.ListListCoset);
    if IndefinitePrint then
        Print("ListCoset=", ListCoset, "\n");
    fi;
    ListCosetRet:=List(ListCoset, x->TheSpaceInv * x * TheSpace);
    return ListCosetRet;
end;


MatrixIntegral_Equivalence_TestFeasibility:=function(GRPrat, EquivRat)
    local TheDenEquiv, TheDenGRP;
    TheDenEquiv:=ReducePrimeMultiplicity(GetDenominatorMatrix(EquivRat));
    TheDenGRP:=GetRationalInvariant(GRPrat);
    if IsInt(TheDenGRP / TheDenEquiv) = false then
        if IndefinitePrint then
            Print("Some prime numbers in the equivalence are not in the group. No equivalence possible");
        fi;
        return false;
    fi;
    return true;
end;


# Find a matrix g in GRPrat such that   g * EquivRat   in   GL(n,Z)
MatrixIntegral_Equivalence:=function(GRPrat, EquivRat)
    local n, LGen, TheSpace, TheSpaceInv, ListGenInt, GRPspace, TheSpaceImg, TheSpaceImgInv, eSpaceEquiv, eMatFinal, eProd;
    if MatrixIntegral_Equivalence_TestFeasibility(GRPrat, EquivRat)=false then
        return fail;
    fi;
    n:=Length(EquivRat);
    LGen:=GeneratorsOfGroup(GRPrat);
    TheSpace:=MatrixIntegral_GetInvariantSpace(n, GRPrat);
    if IndefinitePrint then
        Print("DeterminantMat(TheSpace)=", DeterminantMat(TheSpace), "\n");
    fi;
    # We have TheSpace * g in TheSpace
    # So, in other words TheSpace * g = g_int * TheSpace
    # which gets us TheSpace * g * TheSpaceInv = g_int
    TheSpaceInv:=Inverse(TheSpace);
    ListGenInt:=List(LGen, x->TheSpace * x * TheSpaceInv);
    if First(ListGenInt, x->IsIntegralMat(x)=false)<>fail then
        Error("Some geneatorsmatrix are not integral");
    fi;
    GRPspace:=PersoGroup(ListGenInt, Identity(GRPrat));
    if IsIntegralMat(TheSpaceInv)=false then
        Error("TheSpaceInv is not integral. That should not hapen since TheSpace is built by accretting to Z^n");
    fi;
    # We search for g in GRPrat s.t. g * EquivRat in GL_n(Z).
    # So, we search g in GRPrat s.t. Z^n * g * EquivRat = Z^n
    # Writing g = TheSpaceInv g_int TheSpace we get
    # TheSpaceInv g TheSpace EquivRat = Z^n
    # Or TheSpaceInv g = Inverse(TheSpace * EquivRat)
    TheSpaceImg:=TheSpace * EquivRat;
    TheSpaceImgInv:=Inverse(TheSpaceImg);
    if IsIntegralMat(TheSpaceImgInv)=false then
        return fail;
    fi;
    eSpaceEquiv:=LinearSpace_Equivalence(GRPspace, TheSpaceInv, TheSpaceImgInv);
    if IndefinitePrint then
        Print("Invariant(GRPspace)=", InvariantSpaceOfGroup(n, GRPspace), "\n");
    fi;
    if eSpaceEquiv=fail then
        return fail;
    fi;
    eMatFinal:=TheSpaceInv * eSpaceEquiv * TheSpace;
    eProd:=eMatFinal * EquivRat;
    if IsIntegralMat(eProd)=false then
        Error("The matrix should be integral");
    fi;
    return eProd;
end;


# Find a matrix g in GRPrat such that   EquivRat * g   in   GL(n,Z)
MatrixIntegral_Equivalence_Bis:=function(GRPrat, EquivRat)
    local EquivRatInv, n, TheSol;
    EquivRatInv:=Inverse(EquivRat);
    n:=Length(EquivRat);
    TheSol:=MatrixIntegral_Equivalence(GRPrat, EquivRatInv);
    if TheSol=fail then
        return fail;
    fi;
    # So we have TheSol = g * Inverse(EquivRat) in GL(n,Z)
    # Inverse(TheSol) = EquivRat * g in GL(n,Z)
    return Inverse(TheSol);
end;
