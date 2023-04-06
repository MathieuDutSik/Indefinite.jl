

IsLocallyOKPartialAffineBasis:=function(ThePartial, TheEXT)
  local eVert, B, eVal;
  for eVert in TheEXT
  do
    B:=SolutionMat(ThePartial, eVert);
    if B<>fail then
      for eVal in B
      do
        if IsInt(eVal)=false then
          return false;
        fi;
      od;
    fi;
  od;
  return true;
end;


FindAffineBasisExtension:=function(OldAffineBasis, EXT)
  local eVert, PROV;
  for eVert in EXT
  do
    PROV:=ShallowCopy(OldAffineBasis);
    Add(PROV, eVert);
    if RankMat(PROV)=Length(PROV) then
      if IsLocallyOKPartialAffineBasis(PROV, EXT)=true then
        return PROV;
      fi;
    fi;
  od;
  return false;
end;


Kernel_ExtendToCompleteAffineBasis:=function(EXT, StartingPoint)
  local AffBasis, iRank, TheReply, PersoRank;
  PersoRank:=function(EXT)
    if Length(EXT)=0 then
      return 0;
    fi;
    return RankMat(EXT);
  end;
  AffBasis:=ShallowCopy(StartingPoint);
  for iRank in [PersoRank(AffBasis)+1..PersoRank(EXT)]
  do
    TheReply:=FindAffineBasisExtension(AffBasis, EXT);
    if TheReply=false then
#      Print("afBasis=", AffBasis, "\n");
#      Print("Error, no affine basis found\n");
      return false;
    fi;
    AffBasis:=ShallowCopy(TheReply);
  od;
  return AffBasis;
end;


ExtendToCompleteAffineBasis:=function(EXT, StartingPoint)
  local TheAffBas, len, GRP, nbIter, ePerm, EXTperm;
  TheAffBas:=Kernel_ExtendToCompleteAffineBasis(EXT, StartingPoint);
  if TheAffBas<>false then
    return TheAffBas;
  fi;
  len:=Length(EXT);
  GRP:=SymmetricGroup(len);
  nbIter:=0;
  while(true)
  do
    if nbIter> 10000 then
      Print("Maybe there is no affine basis after all");
      return false;
    fi;
    ePerm:=Random(GRP);
    EXTperm:=Permuted(EXT, ePerm);
    TheAffBas:=Kernel_ExtendToCompleteAffineBasis(EXTperm, StartingPoint);
    if TheAffBas<>false then
      return TheAffBas;
    fi;
    nbIter:=nbIter+1;
  od;
end;
