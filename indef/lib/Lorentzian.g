

# This version works only if Subspace1 is a codimension 1 subspace.
LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1:=function(G1, Subspace1, G2, Subspace2)
    local dim, Compl1, eVect1, eNorm, eProd1, eProd2, Vscal, V0, NSP, V1, eNorm_V1, scal0, scal1, tVal, eVect2, Trans1, Trans2, eEquiv, eInv, G1_tr;
    dim:=Length(G1);
    Compl1:=SubspaceCompletionRational(Subspace1, dim);
    if Length(Compl1)<>1 then
        Error("Compl1 should be of length 1");
    fi;
    eVect1:=Compl1[1];
    eNorm:=eVect1 * G1 * eVect1;
    eProd1 := Subspace1 * G1;
    eProd2 := Subspace2 * G2;
    Vscal:=List(eProd1, x->x*eVect1);
    V0:=SolutionMat(TransposedMat(eProd2), Vscal);
    if V0=fail then
        Error("The solutionning failed");
    fi;
    NSP:=NullspaceMat(TransposedMat(eProd2));
    if Length(NSP)<>1 then
        Error("NSP should have length 1");
    fi;
    V1:=NSP[1];
    eNorm_V1 := V1 * G2 * V1;
    if eNorm_V1<>0 then
        Error("V1 should be an isotrop vector");
    fi;
    scal0 := eNorm - V0 * G2 * V0;
    scal1 := 2 * (V0 * G2 * V1);
    if scal1=0 then
        Error("scal1 should be non-zero");
    fi;
    tVal:=scal0 / scal1;
    eVect2:=V0 + tVal * V1;
    Trans1:=Concatenation(Subspace1, [eVect1]);
    Trans2:=Concatenation(Subspace2, [eVect2]);
    eEquiv:=Inverse(Trans1) * Trans2;
    eInv:=Inverse(eEquiv);
    G1_tr:=eInv * G1 * TransposedMat(eInv);
    if G1_tr<>G2 then
        Error("G1 was not mapped to G2");
    fi;
    if Subspace1*eEquiv <> Subspace2 then
        Error("Subspace1 is not mapped to Subspace2");
    fi;
    return eEquiv;
end;


# Resolution of B = X A + A^T X^T
# b_{ij} = sum_k x_{ik} a_{kj} + a_{ki} x_{jk}
SpecialEquationSolving:=function(Amat, Bmat)
    local dim, TheMat, Bvect, f, i, j, u, k, v1, v2, f_getmat, eSol_vect, eSol_mat, BasisKernel, NSP, eVect, eMat;
    dim:=Length(Bmat);
    TheMat:=NullMat(dim * dim, dim * dim);
    Bvect:=ListWithIdenticalEntries(dim * dim, 0);
    f:=function(i, j)
        return i + dim * (j-1);
    end;
    for i in [1..dim]
    do
        for j in [1..dim]
        do
            u:=f(i,j);
            Bvect[u] := Bmat[i][j];
            for k in [1..dim]
            do
                # contribudion x_{ik} a_{kj}
                v1:=f(i, k);
                TheMat[v1][u] := TheMat[v1][u] + Amat[k][j];
                # contribution a_{ki} x_{jk}
                v2:=f(j, k);
                TheMat[v2][u] := TheMat[v2][u] + Amat[k][i];
            od;
        od;
    od;
    f_getmat:=function(eVect)
        local eMat, i, j;
        eMat:=NullMat(dim, dim);
        for i in [1..dim]
        do
            for j in [1..dim]
            do
                eMat[i][j] := eVect[f(i,j)];
            od;
        od;
        return eMat;
    end;
    eSol_vect:=SolutionMat(TheMat, Bvect);
    eSol_mat:=f_getmat(eSol_vect);
    if Bmat<>eSol_mat * Amat + TransposedMat(Amat) * TransposedMat(eSol_mat) then
        Error("Failed to find the correct solution 1");
    fi;
    BasisKernel:=[];
    NSP:=NullspaceMat(TheMat);
    for eVect in NSP
    do
        eMat:=f_getmat(eVect);
        if NullMat(dim,dim)<>eMat * Amat + TransposedMat(Amat) * TransposedMat(eMat) then
            Error("Failed to find the correct solution 2");
        fi;
        Add(BasisKernel, eMat);
    od;
    return rec(BasisKernel:=BasisKernel, eSol_mat:=eSol_mat);
end;


# Subspace1 and Subspace2 are in the orthogonal of a K-space of isotropic vectors.
# Subspace1 is a list of vectors and Subspace2 another list of vectors.
# We build a full dimensional orthogonal mapping that maps 1 on 2.
# and that is an orthogonal transformation mapping G1 to G2.
#
# The returned extension is not unique if Subspace1 / Subspace2 have codimension > 1.
LORENTZ_ExtendOrthogonalIsotropicIsomorphism:=function(G1, Subspace1, G2, Subspace2)
    local dim, rnk, dimSpace, dimCompl, TheCompl1, TheCompl2, eVect1, eProd1, eProd2, Vscal, NSP1, NSP2, ProdMat1, ProdMat2, eVect2, Trans1, Trans1Inv, Trans2, eEquiv, eEquivInv, MatScal1, MatScal2, MatScalCand2, G1_tr, ListVectCand2, DiffScal, LambdaMat, LScal1, LScal2, eVectCand2, SMat1, SMat2, denomSol, test_transformation, check_transformation, get_one_transformation, get_kernel_generating_set;
    dim:=Length(G1);
    rnk:=RankMat(Subspace1);
#    Print("------------------------------\n");
#    Print("Subspace1:=", Subspace1, ";\n");
#    Print("Subspace2:=", Subspace2, ";\n");
    if Length(Subspace1)<>rnk then
        Error("Inconsistent input: We should have Subspace1 of full rank");
    fi;
    dimSpace:=rnk;
    dimCompl:=dim - dimSpace;
#    Print("dim:=", dim, ", dimSpace:=", dimSpace, ", dimCompl:=", dimCompl, ";\n");
    # Checking the input
    eProd1:=Subspace1 * G1;
    eProd2:=Subspace2 * G2;
    NSP1:=NullspaceMat(TransposedMat(eProd1));
    ProdMat1:=NSP1 * G1 * TransposedMat(NSP1);
    if ProdMat1<>NullMat(dimCompl,dimCompl) then
        Error("Inconsistent input: ProdMat1 should be equal to zero");
    fi;
    NSP2:=NullspaceMat(TransposedMat(eProd2));
    ProdMat2:=NSP2 * G2 * TransposedMat(NSP2);
    if ProdMat2<>NullMat(dimCompl,dimCompl) then
        Error("Inconsistent input: ProdMat2 should be equal to zero");
    fi;
    SMat1:=Subspace1 * G1 * TransposedMat(Subspace1);
    SMat2:=Subspace2 * G2 * TransposedMat(Subspace2);
    if SMat1<>SMat2 then
        Error("Inconsistent input: SMat1 should be equal to SMat2");
    fi;
    # Now doing the computation side of the job
    TheCompl1:=SubspaceCompletionRational(Subspace1, dim);
    Trans1:=Concatenation(Subspace1, TheCompl1);
    Trans1Inv:=Inverse(Trans1);
    MatScal1:=TheCompl1 * G1 * TransposedMat(TheCompl1);
    ListVectCand2:=[];
    for eVect1 in TheCompl1
    do
        Vscal:=List(eProd1, x->x*eVect1);
        eVectCand2:=SolutionMat(TransposedMat(eProd2), Vscal);
        if eVectCand2 = fail then
            Error("The solutionning failed");
        fi;
        LScal1:=List(Subspace1, x->x * G1 * eVect1);
        LScal2:=List(Subspace2, x->x * G2 * eVectCand2);
        if LScal1<>LScal2 then
            Error("Inconsistency for LScal1 = LScal2");
        fi;
        Add(ListVectCand2, eVectCand2);
    od;
    # The solutions are written as
    # eVect2 = eVectCand2 + c_vect * NSP2
    # Putting together this gets ListVect2 = TheCompl2 = ListVectCand2 + c_Mat * NSP2
    # MatScal2 = (ListVectCand2 + c_Mat * NSP2) * G2 * (NSP2^T * c_Mat^T + ListVectCand2^T)
    #          = ListVectCand2 * G2 * ListVectCand2^T + c_Mat * NSP2 * G2 * ListVectCand2^T + ListVectCand2 * G2 * NSP2^T * c_Mat^T
    # c_vect is a vector of length dimCompl which gets into
    # Put all together, this gets us c_Mat a matrix of size (dimCompl, dimCompl)
    # The equation that we get is thus of the form B = X A + A^T X^T
    # The equation is underdefined. This is because B is symmetric and so we have n(n+1)/2 equations
    # for n^2 unknowns.
    # The space NSP2 is uniquely defined as the set of isotropic vectors in the space. It is defined
    # as a Kernel.
    MatScalCand2:=ListVectCand2 * G2 * TransposedMat(ListVectCand2);
    LambdaMat := NSP2 * G2 * TransposedMat(ListVectCand2);
    denomSol:=GetDenominatorQuotientSolution(LambdaMat);
    DiffScal:=MatScal1 - MatScalCand2;
    test_transformation:=function(eEq)
        local eEqInv, G1_tr;
        eEqInv:=Inverse(eEq);
        G1_tr:=eEqInv * G1 * TransposedMat(eEqInv);
        if G1_tr<>G2 then
#            Print("G1 was not mapped to G2\n");
            return false;
        fi;
        if Subspace1*eEq <> Subspace2 then
#            Print("Subspace1 is not mapped to Subspace2\n");
            return false;
        fi;
        return true;
    end;
    check_transformation:=function(eEq)
        if test_transformation(eEq)=false then
            Error("Something wrong happened");
        fi;
    end;
    get_one_transformation:=function()
        local TheRec, f_get_equiv, eEquiv0, ListEquiv_terms, eEquiv0_v, ListEquiv_terms_v, TheSol, RetSol;
        TheRec:=SpecialEquationSolving(LambdaMat, DiffScal);
        f_get_equiv:=function(eSol)
            local TheCompl2, Trans2;
            TheCompl2:=ListVectCand2 + eSol * NSP2;
            Trans2:=Concatenation(Subspace2, TheCompl2);
            return Trans1Inv * Trans2;
        end;
        eEquiv0:=f_get_equiv(TheRec.eSol_mat);
        ListEquiv_terms:=List(TheRec.BasisKernel, x->Trans1Inv * Concatenation(NullMat(rnk, dim), x * NSP2));
        # The matrix is expressed as eEquiv0 + alpha1 ListEquiv_terms[1] + ..... + alphaN ListEquiv_terms[N]
        RetSol:=EliminateSuperfluousPrimeDenominators_Matrix(eEquiv0, ListEquiv_terms);
        check_transformation(RetSol);
        return RetSol;
    end;
    get_kernel_generating_set:=function(d)
        local TheRec, ListEquiv_terms1, ListEquiv_terms2, ListEquiv_terms3, ListEquiv_terms4, ListEquiv_terms5, eGen;
#        Print("Beginning of get_kernel_generating_set\n");
        if G1<>G2 or Subspace1<>Subspace2 then
            Error("We should have G1=G2 and Subspace1=Subspace2 in order for kernel to make sense");
        fi;
        TheRec:=SpecialEquationSolving(LambdaMat, DiffScal);
        ListEquiv_terms1:=List(TheRec.BasisKernel, x->Trans1Inv * Concatenation(NullMat(rnk, dim), x * NSP2));
        ListEquiv_terms2:=List(ListEquiv_terms1, MatrixToVector);
        ListEquiv_terms3:=IntegralSpaceSaturation(ListEquiv_terms2);
        ListEquiv_terms4:=List(ListEquiv_terms3, VectorToMatrix);
        ListEquiv_terms5:=List(ListEquiv_terms4, x->IdentityMat(dim) + x / d);
        for eGen in ListEquiv_terms5
        do
            check_transformation(eGen);
        od;
        return PersoGroupMatrix(ListEquiv_terms5, dim);
    end;
    return rec(denomSol:=denomSol,
               check_transformation:=check_transformation,
               get_one_transformation:=get_one_transformation,
               get_kernel_generating_set:=get_kernel_generating_set);
end;


#
# this is a set of functions related to Lorentzian lattice computation
# in one way or another.
LORENTZ_IsLorentzian:=function(LorMat)
  local eRec;
  eRec:=DiagonalizeSymmetricMatrix(LorMat);
  if eRec.nbPlus<>1 then
#    Print("Cannot be lorentzian. Plus space should be one-dimensional\n");
    return false;
  fi;
  if eRec.nbZero<>0 then
#    Print("Cannot be lorentzian. Should be nondegenerate\n");
    return false;
  fi;
  return true;
end;


LORENTZ_AssembleDiagBlocks:=function(ListMat)
  local nbMat, ListDim, TotDim, RetMat, iMat, TheShift, jMat, eDim, i, j;
  nbMat:=Length(ListMat);
  ListDim:=List(ListMat, Length);
  TotDim:=Sum(ListDim);
  RetMat:=NullMat(TotDim, TotDim);
  for iMat in [1..nbMat]
  do
    TheShift:=0;
    for jMat in [1..iMat-1]
    do
      TheShift:=TheShift + ListDim[jMat];
    od;
    eDim:=ListDim[iMat];
    for i in [1..eDim]
    do
      for j in [1..eDim]
      do
        RetMat[i+TheShift][j+TheShift]:=ListMat[iMat][i][j];
      od;
    od;
  od;
  return RetMat;
end;


#
#
# We check that all the vectors are in a plane
LORENTZ_CheckCorrectnessVectorFamily:=function(ListVect)
  local dim, ListVectExt, rnk;
  dim:=Length(ListVect[1]);
  ListVectExt:=List(ListVect, x->Concatenation([1], x));
  rnk:=RankMat(ListVectExt);
  if rnk<> dim then
    Print("We have dim=", dim, " rnk=", rnk, " they should be equal\n");
    Error("Please correct this problem");
  fi;
end;


#
# We enumerate all the vectors fVect such that
# fVect * LorMat * fVect >= 0 (or =0)     and      0 < fVect * LorMat * eVect <= MaxScal
#
# TheOption can be "isotrop" or "total"
# isotrop means all the vectors that are isotrop
#
# OnlyShortest is true if we return only the vectors of minimal values fVect * LorMat * eVect
LORENTZ_FindPositiveVectors:=function(LorMat, eVect, MaxScal, TheOption, OnlyShortest)
  local eNorm, Ubasis, GramMat, TheRec, eValMax, TotalListSol, eVal, eBasSol, alpha, eTrans, eSol, eSquareDist, ListSol, eSolA, eSolB, eSolC, fScal, fNorm;
#  Print("Beginning of LORENTZ_FindPositiveVectors\n");
  eNorm:=eVect*LorMat*eVect;
  if MaxScal<=0 then
    Print("MaxScal=", MaxScal, "\n");
    Error("We should have MaxScal > 0");
  fi;
  if TheOption<>"isotrop" and TheOption<>"total" then
    Print("TheOption=", TheOption, "\n");
    Error("AllowedValues for TheOption are isotrop or total");
  fi;
  if eNorm <= 0 then
    Print("eNorm=", eNorm, "\n");
    Error("Wrong norm of vector, will not work 1");
  fi;
  if IsIntegralVector(eVect)=false then
    Print("eVect=", eVect, "\n");
    Error("Vector should be integral.");
  fi;
  Ubasis:=NullspaceIntMat(TransposedMat([eVect*LorMat]));
  GramMat:=-Ubasis*LorMat*TransposedMat(Ubasis);
  if IsPositiveDefiniteSymmetricMatrix(GramMat)=false then
    Error("deep error");
  fi;
  TheRec:=GcdVector(eVect*LorMat);
  if TheRec.TheGcd<0 then
    Error("deep error 2");
  fi;
  eValMax:=LowerInteger(MaxScal/TheRec.TheGcd);
#  Print("MaxScal=", MaxScal, " TheRec.TheGcd=", TheRec.TheGcd, " eValMax=", eValMax, "\n");
  TotalListSol:=[];
  for eVal in [1..eValMax]
  do
    eBasSol:=eVal*TheRec.ListCoef; # a solution of eVect * LorMat * fVect = TheGcd * eVal
    alpha:=eVal*TheRec.TheGcd/eNorm;
    eTrans:=eBasSol - alpha*eVect;
    eSol:=SolutionMat(Ubasis, eTrans);
    if eSol=fail then
      Error("Please debug from here Lorentian 1");
    fi;
    eSquareDist:=alpha*alpha*eNorm;
    ListSol:=ClosestAtDistanceVallentinProgram(GramMat, -eSol, eSquareDist);
    for eSolA in ListSol
    do
      eSolB:=eSolA*Ubasis;
      eSolC:=eBasSol+eSolB;
      fScal:=eSolC*LorMat*eVect;
      fNorm:=eSolC*LorMat*eSolC;
      if fNorm < 0 or fScal < 0 or fScal > MaxScal then
        Error("Error in norm or scalar product");
      fi;
      if IsIntegralVector(eSolC)=false then
        Error("Non integral eSolC, debug ....");
      fi;
      if TheOption="total" then
        Add(TotalListSol, eSolC);
      else
        if fNorm=0 then
          Add(TotalListSol, eSolC);
        fi;
      fi;
    od;
    if OnlyShortest and Length(TotalListSol)>0 then
      return TotalListSol;
    fi;
  od;
  return TotalListSol;
end;


LORENTZ_SearchInitialVectorOpt:=function(LorMat, PosVect, TheOption, OnlyShortest)
  local n, MaxScal, ListPosVect, ListIsotrop, ListScal, eScal, eSet;
  n:=Length(LorMat);
  MaxScal:=PosVect*LorMat*PosVect;
#  Print("Beginning of LORENTZ_SearchInitialVector\n");
#  Print("LorMat=\n");
#  PrintArray(LorMat);
#  Print("PosVect=", PosVect, " MaxScal=", MaxScal, "\n");
  while(true)
  do
    ListPosVect:=LORENTZ_FindPositiveVectors(LorMat, PosVect, MaxScal, TheOption, OnlyShortest);
#    Print("LORENTZ_SearchInitialVector MaxScal=", MaxScal, " |ListPosVect|=", Length(ListPosVect), "\n");
    if Length(ListPosVect)>0 then
      ListScal:=List(ListPosVect, x->x*LorMat*PosVect);
      eScal:=Minimum(Difference(Set(ListScal), [0]));
      eSet:=Filtered([1..Length(ListScal)], x->ListScal[x]=eScal);
      return ListPosVect{eSet};
    fi;
    MaxScal:=MaxScal+1;
  od;
end;


LORENTZ_SearchInitialVector:=function(LorMat, PosVect, TheOption)
    local OnlyShortest;
    OnlyShortest:=false;
    return LORENTZ_SearchInitialVectorOpt(LorMat, PosVect, TheOption, OnlyShortest);
end;


LORENTZ_GetDefiningIneq:=function(LorMat, ListVect)
  local nbVect, ListB, eSol;
  nbVect:=Length(ListVect);
  ListB:=ListWithIdenticalEntries(nbVect,1);
  eSol:=SolutionMat(TransposedMat(ListVect), ListB);
  return RemoveFraction(eSol*Inverse(LorMat));
end;


LORENTZ_GetDeterminingVectFamily:=function(LorMat, eFamEXT)
  local eVectDef, ListScal, MaxScal, ListVect, TheDet, AffBas, TheSel, ListCollScal, LColl, LSize, LType, nbCase, iCase, ePerm, LTypeS, ListVectB, CurrDet, eListSel, TestVect, NewDet, OnlyShortest;
#  Print("Running LORENTZ_GetDeterminingVectFamily |eFamEXT|=", Length(eFamEXT), "\n");
  eVectDef:=LORENTZ_GetDefiningIneq(LorMat, eFamEXT);
  ListScal:=List(eFamEXT, x->x*LorMat*eVectDef);
  if Length(Set(ListScal))<>1 then
    Error("We are wrong");
  fi;
  MaxScal:=ListScal[1];
  TheDet:=AbsInt(DeterminantMat(BaseIntMat(eFamEXT)));
#  Print("MaxScal=", MaxScal, " det=", TheDet, "\n");
  if TheDet=1 then
    return eFamEXT;
  fi;
  while(true)
  do
    OnlyShortest:=false;
    ListVect:=LORENTZ_FindPositiveVectors(LorMat, eVectDef, MaxScal, "total", OnlyShortest);
    TheDet:=AbsInt(DeterminantMat(BaseIntMat(ListVect)));
    if TheDet=1 then
      AffBas:=GetZbasis(eFamEXT);
      TheSel:=Filtered(ListVect, x->SolutionIntMat(AffBas, x)=fail);
      ListCollScal:=List(TheSel, x->Collected(List(eFamEXT, y->y*LorMat*x)));
      LColl:=Collected(ListCollScal);
      LSize:=List(LColl, x->x[2]);
      LType:=List(LColl, x->x[1]);
      nbCase:=Length(LColl);
      ePerm:=SortingPerm(LSize);
      LTypeS:=Permuted(LType, ePerm);
      ListVectB:=ShallowCopy(eFamEXT);
      CurrDet:=AbsInt(DeterminantMat(BaseIntMat(eFamEXT)));
      for iCase in [1..nbCase]
      do
        eListSel:=Filtered([1..Length(TheSel)], x->ListCollScal[x]=LTypeS[iCase]);
        TestVect:=Concatenation(ListVectB, TheSel{eListSel});
        NewDet:=AbsInt(DeterminantMat(BaseIntMat(TestVect)));
        if NewDet<CurrDet then
          ListVectB:=ShallowCopy(TestVect);
          CurrDet:=NewDet;
          if CurrDet=1 then
            return ListVectB;
          fi;
        fi;
      od;
    fi;
    MaxScal:=MaxScal+1;
  od;
end;


LORENTZ_ComputeFundamentalStabInfo:=function(LorMat, eFamEXT, GRPint)
  local ListPermGen, ListMatrGenB, eGen, eMatr, GRPintMatr, PhiPermMat;
  ListPermGen:=SmallGeneratingSet(GRPint);
  ListMatrGenB:=[];
  for eGen in ListPermGen
  do
    eMatr:=FindTransformation(eFamEXT, eFamEXT, eGen);
    if IsIntegralMat(eMatr)=false then
      Error("Bug: Non integral matrix");
    fi;
    if eMatr*LorMat*TransposedMat(eMatr)<>LorMat then
      Error("Bug: eMatr does not preserve LorMat");
    fi;
    Add(ListMatrGenB, eMatr);
  od;
  GRPintMatr:=Group(ListMatrGenB);
  PhiPermMat:=GroupHomomorphismByImagesNC(GRPint, GRPintMatr, ListPermGen, ListMatrGenB);
  return rec(GRP_int:=GRPint, GRPintMatr:=GRPintMatr, PhiPermMat:=PhiPermMat);
end;


LORENTZ_ComputeStabilizer:=function(LorMat, eFamEXT)
  local GRPrat, ListMatrGen, test, ListVect, GRPtot, ListPermGen, ListMatrGenB, eGen, eMatr, eList, GRPint, GRPintMatr, eVect, testB, ListPGen, PhiPermMat;
#  Print("|eFamEXT|=", Length(eFamEXT), "\n");
  GRPrat:=LinPolytope_Automorphism_AddMat(eFamEXT, [LorMat]);
#  Print("|GRPrat|=", Order(GRPrat), "\n");
  ListPGen:=GeneratorsOfGroup(GRPrat);
#  Print("Before ListMatrGen construction\n");
  ListMatrGen:=List(ListPGen, x->FindTransformation(eFamEXT, eFamEXT, x));
#  Print("After ListMatrGen construction\n");
  test:=First(ListMatrGen, x->IsIntegralMat(x)=false);
  if test=fail then
    GRPintMatr:=Group(ListMatrGen);
    PhiPermMat:=GroupHomomorphismByImagesNC(GRPrat, GRPintMatr, ListPGen, ListMatrGen);
    return rec(GRP_rat:=GRPrat, GRP_int:=GRPrat, GRPintMatr:=GRPintMatr, PhiPermMat:=PhiPermMat);
  fi;
  GRPint:=KernelLinPolytopeIntegral_Automorphism_Subspaces(eFamEXT, GRPrat).GRPperm;
  return LORENTZ_ComputeFundamentalStabInfo(LorMat, eFamEXT, GRPint);
end;


LORENTZ_GetAnsatzGraphInformation:=function(LorMat, eFamEXT, Qmat, HeuristicScal)
  local nbVert, ThePartition, TheListAdjacency, i, eList, j, eScal, ScalarMat, DistMat, korig, pos, LineScalar, iVert, jVert, SetV, n, SecVal;
  nbVert:=Length(eFamEXT);
#  Print("Computing ScalarMat for nbVert=", nbVert, "\n");
  ScalarMat:=[];
  for i in [1..nbVert]
  do
    LineScalar:=[];
    for j in [1..nbVert]
    do
      eScal:=[];
      if HeuristicScal.UseLorMat then
        Add(eScal, eFamEXT[i] * LorMat * eFamEXT[j]);
      fi;
      if HeuristicScal.UseQmat then
        Add(eScal, eFamEXT[i] * Qmat * eFamEXT[j]);
      fi;
      if HeuristicScal.UseAllValues then
        Add(LineScalar, eScal);
      else
        pos:=Position(HeuristicScal.ListAllowedValues, eScal);
        Add(LineScalar, eScal);
      fi;
    od;
    Add(ScalarMat, LineScalar);
  od;
#  Print("We have ScalarMat\n");
  if HeuristicScal.UseDiagonal then
#    Print("Using the diagonal\n");
    DistMat:=MappedScalarMatrixDistanceMatrix(ScalarMat);
    SetV:=__SetValue(DistMat);
    korig:=Length(SetV);
    n:=Length(DistMat);
    TheListAdjacency:=Method4modelEdgeColoredGraph(DistMat, SetV);
    ThePartition:=__Method4Partition(korig, n);
    return rec(ThePartition:=ThePartition, ListAdjacency:=TheListAdjacency);
  else
    SetV:=__SetValue(ScalarMat);
#    Print("Not using the diagonal |SetV|=", Length(SetV), "\n");
    if Length(SetV)=2 then
      SecVal:=SetV[2];
      ThePartition:=[[1..nbVert]];
      TheListAdjacency:=[];
      for iVert in [1..nbVert]
      do
        eList:=[];
        for jVert in [1..nbVert]
        do
          if iVert<>jVert then
            if ScalarMat[iVert][jVert]=SecVal then
              Add(eList, jVert);
            fi;
          fi;
        od;
        Add(TheListAdjacency, eList);
      od;
      return rec(ThePartition:=ThePartition, ListAdjacency:=TheListAdjacency);
    else
      SetV:=__SetValue(ScalarMat);
      korig:=Length(SetV);
      n:=Length(ScalarMat);
      TheListAdjacency:=Method4modelEdgeColoredGraph(ScalarMat, SetV);
      ThePartition:=__Method4Partition(korig, n);
      return rec(ThePartition:=ThePartition, ListAdjacency:=TheListAdjacency);
    fi;
  fi;
end;


LORENTZ_GetAnsatzSubpolytope:=function(LorMat, eFamEXT, Qmat, HeuristicSub)
  local nbVert, TheSub, iVert, eScal;
  nbVert:=Length(eFamEXT);
  TheSub:=[];
  for iVert in [1..nbVert]
  do
    eScal:=[];
    if HeuristicSub.UseLorMat then
      Add(eScal, eFamEXT[iVert] * LorMat * eFamEXT[iVert]);
    fi;
    if HeuristicSub.UseQmat then
      Add(eScal, eFamEXT[iVert] * Qmat * eFamEXT[iVert]);
    fi;
    if HeuristicSub.UseAllValues then
      Add(TheSub, iVert);
    else
      if Position(HeuristicSub.ListAllowedValues, eScal)<>fail then
        Add(TheSub, iVert);
      fi;
    fi;
  od;
  return TheSub;
end;


LORENTZ_TestEquivalence_CheckAndReturn:=function(LorMat1, LorMat2, eFamEXT1, eFamEXT2, eMatrB)
  if IsIntegralMat(eMatrB)=false then
    Error("Bug: he matrix should be integral");
  fi;
  if Inverse(eMatrB)*LorMat1*TransposedMat(Inverse(eMatrB))<>LorMat2 then
    Error("Bug: the matrix does not map the Lorentzian quadratic form");
  fi;
  if Set(eFamEXT1*eMatrB)<>Set(eFamEXT2) then
    Error("Bug: the matrix does not map the vector configurations");
  fi;
  return eMatrB;
end;


LORENTZ_TestEquivalence_General:=function(LorMat1, LorMat2, eFamEXT1, eFamEXT2)
  local eEquiv, eEquivB, eMatr, TheStrategy, ListVect1, ListVect2, eMatrB, GRP2;
#  Print("Before LinPolytope_Isomorphism_AddMat |eFamEXT1|=", Length(eFamEXT1), " |eFamEXT2|=", Length(eFamEXT2), "\n");
  eEquiv:=LinPolytope_Isomorphism_AddMat(eFamEXT1, eFamEXT2, [LorMat1], [LorMat2]);
#  Print("After LinPolytope_Isomorphism_AddMat\n");
  if eEquiv=false then
    return false;
  fi;
  eMatr:=FindTransformation(eFamEXT1, eFamEXT2, eEquiv);
  if IsIntegralMat(eMatr) and Inverse(eMatr)*LorMat1*TransposedMat(Inverse(eMatr))=LorMat2 then
    return eMatr;
  fi;
  if AbsInt(DeterminantMat(BaseIntMat(eFamEXT1)))=1 then
    Error("It should have ended at this stage");
  fi;
#  Print("Before LinPolytope_Automorphism_AddMat\n");
  GRP2:=LinPolytope_Automorphism_AddMat(eFamEXT2, [LorMat2]);
#  Print("After LinPolytope_Automorphism_AddMat\n");
  eMatrB:=KernelLinPolytopeIntegral_Isomorphism_Subspaces(eFamEXT1, eFamEXT2, GRP2, eEquiv);
  if eMatrB=false then
#    Print("Subspaces algo returns false\n");
    return false;
  fi;
#  Print("We have eMatrB\n");
  return LORENTZ_TestEquivalence_CheckAndReturn(LorMat1, LorMat2, eFamEXT1, eFamEXT2, eMatrB);
end;


LORENTZ_TestEquivalence:=function(LorMat, eFamEXT1, eFamEXT2)
  return LORENTZ_TestEquivalence_General(LorMat, LorMat, eFamEXT1, eFamEXT2);
end;


GetRationalIsotropyVectors:=function(LorMat22)
    local a, b, c, DeltaB, TheSqrt, root1, root2, v1, v2;
    # The matrix is written as
    # | a b |
    # | b c |
    # The quadratic form is a x^2 + 2b x y + c y^2
    # We have Delta = 4 (b^2 - a c)
    # So { a (x / y)^2 + 2b (x/y) + c } y^2
    # This gets us
    # x1 / y1 = (-2 b + sqrt(Delta) ) / (2a)
    #         = (-b + sqrt(DeltaB)) / a
    # x2 / y2 = (-2 b - sqrt(Delta) ) / (2a)
    #         = (-b - sqrt(DeltaB) ) / a
    a:=LorMat22[1][1];
    b:=LorMat22[1][2];
    c:=LorMat22[2][2];
    DeltaB:=b*b - a*c;
    if DeltaB < 0 then
        Error("If the discriminant is negative, then we cannot have signature (1,1)");
    fi;
    TheSqrt:=RatSqrt(DeltaB);
    if TheSqrt = fail then
        return [];
    fi;
    if a<>0 then
        root1:=(-b + TheSqrt) / a;
        root2:=(-b - TheSqrt) / a;
        v1:=[ root1, 1];
        v2:=[ root2, 1];
        return [v1, v2];
    fi;
    v1:=[1,0];
    # Remaining equation is 2 b x + c y = 0
    v2:=[-c, 2*b];
    return [v1, v2];
end;


# We compute an upper bound TheMult on the possible values of x such that
# eNSP = eNSPbas + x * eNSPdir is admissible as a facet vector of the
#
# We also check for the isotropy situation that could happen and report it
# on the output. That is if the best x_upp gives eNSP[1] = 0 and eNSP{[2..n+1]}
# is isotropic than in the following call to the LORENTZ_FindPositiveVectors
# the following happens:
# -- MaxScal = 0 (because eNSP[1] = 0)
# -- eVect * LorMat * fVect = 0 (because 0 <= val <= 0)
# -- Thus is fVect is in the orthogonal of an isotropic vector and is of norm >= 0.
# -- By the geometry this gets us to fVect a multiple of eVect
# That scenario is not acceptable for finding perfect domain.
GetUpperBound:=function(LorMat, eNSPbas, eNSPdir)
    local n, eCstBas, eCstDir, eBas, eDir, eVectBas, eVectDir, eNorm, iShift, eV, eVect, ListUpperBound, TheBasis, LorMat22, fact, ListIso, eIso, BestUpper, UpperBound_constant, UpperBound_isotropic, ListUpperBound_Iso;
    n:=Length(LorMat);
    eCstBas:=eNSPbas[1];
    eCstDir:=eNSPdir[1];
    eBas:=eNSPbas{[2..n+1]};
    eDir:=eNSPdir{[2..n+1]};
#    Print("eCstBas=", eCstBas, " eCstDir=", eCstDir, "\n");
    # For an acceptable vector w of length n+1 in return we must have w[1] < 0.
    # Since we have w = u + TheMult * v we have a potential upper bound
    # on TheMult, but only if v[1] > 0
    ListUpperBound:=[];
    UpperBound_constant:=fail;
    if eCstDir>0 then
        UpperBound_constant:=-eCstBas / eCstDir;
        if UpperBound_constant <= 0 then
            Error("The upper bound from constant is not correct");
        fi;
        Add(ListUpperBound, UpperBound_constant);
    fi;
    #
    # Get Raw upper bound
    #
    iShift:=1;
    while(true)
    do
        eV:=eBas + iShift * eDir;
        eVect:=eV * Inverse(LorMat);
        if eVect * LorMat * eVect < 0 then
            Add(ListUpperBound, iShift);
            break;
        fi;
        iShift:=iShift * 2;
    od;
    #
    # More subttle upper bound coming from isotropy computation
    #
    eVectBas:=eBas * Inverse(LorMat);
    eVectDir:=eDir * Inverse(LorMat);
    TheBasis:=[eVectBas, eVectDir];
    LorMat22:=TheBasis * LorMat * TransposedMat(TheBasis);
    ListIso:=GetRationalIsotropyVectors(LorMat22);
    UpperBound_isotropic:=fail;
    ListUpperBound_Iso:=[];
    for eIso in ListIso
    do
        if eIso[1] <> 0 then
            fact:=eIso[2] / eIso[1];
            eVect:=eVectBas + fact * eVectDir;
            if eVect * LorMat * eVect <> 0 then
                Error("eVect should be isotropic");
            fi;
            if fact > 0 then
                Add(ListUpperBound_Iso, fact);
            fi;
        fi;
    od;
    if Length(ListUpperBound_Iso) > 0 then
        UpperBound_isotropic:=Minimum(ListUpperBound_Iso);
        Add(ListUpperBound, UpperBound_isotropic);
    fi;
    if UpperBound_constant<>fail and UpperBound_isotropic<>fail then
        if UpperBound_constant = UpperBound_isotropic then
            return fail;
        fi;
    fi;
    # Need to see if better upper bounds are possible, but this is a secondary question
    BestUpper:=Minimum(ListUpperBound);
    eV:=eBas + BestUpper * eDir;
    eVect:=eV * Inverse(LorMat);
    eNorm:=eVect * LorMat * eVect;
    if eNorm > 0 then
        Error("We should have eNorm <= 0");
    fi;
    return BestUpper;
end;


# Given a Critical set a vector eNSPbas and a direction eNSPdir
# we are looking for a lambda > 0 such that eNSP = eNSPbas + lambda * eNSPdir
# such that the list of vectors satisfying eNSP * v = 0 defines a bigger
# set than CritSet.
#
# eNSPbas must be of positive norm. eNSPdir must be of negative norm.
#
LORENTZ_Kernel_Flipping:=function(LorMat, CritSet, eNSPbas, eNSPdir, TheOption)
  local n, EXT, NSP, eVert, eEXT, eNSP, NSPb, eNSPb, eNSPtest, eVectTest, MaxScal, ListTotal, EXTtest, ListIsoTest, NSPtest, eVectB, eNormTest, eVect, OnlyShortest, eDen, aShift, TheLowerBound, TheUpperBound, uVal, vVal, TheMidVal, n_iter, eVectBas, eVectDir, eNormBas, eNormDir, eVectBasDir, eNormBasDir, GetMidVal;
#  Print("Beginning of LORENTZ_Kernel_Flipping\n");
#  Print("LorMat=\n");
#  PrintArray(LorMat);
#  Print("CritSet=");
#  PrintArray(CritSet);
#  Print("eNSPbas=", eNSPbas, "  eNSPdir=", eNSPdir, "\n");
#  Print("TheOption=", TheOption, "\n");
  if RankMat([eNSPbas,eNSPdir])<>2 then
      Error("The vector eNSPbas and eNSPdir should be linearly independent");
  fi;
  if First(CritSet, x->Concatenation([1],x)*eNSPbas<>0)<>fail then
      Error("eNSPbas should have scalar product 0 with all entries in CritSet");
  fi;
  if First(CritSet, x->Concatenation([1],x)*eNSPdir<>0)<>fail then
      Error("eNSPdir should have scalar product 0 with all entries in CritSet");
  fi;
  if TheOption="isotrop" then
    if First(CritSet, x->x*LorMat*x<>0)<>fail then
      Print("CritSet=", CritSet, "\n");
      Error("CritSet contains some non-isotrop vectors");
    fi;
  fi;
  n:=Length(LorMat);
  OnlyShortest:=true;
  TheLowerBound:=0;
  TheUpperBound:=GetUpperBound(LorMat, eNSPbas, eNSPdir);
#  Print("At the beginning TheUpperBound=", TheUpperBound, "\n");
  n_iter:=0;
  eVectBas:=RemoveFraction(eNSPbas{[2..n+1]}*Inverse(LorMat));
  eVectDir:=RemoveFraction(eNSPdir{[2..n+1]}*Inverse(LorMat));
  eNormBas:=eVectBas * LorMat * eVectBas;
  eNormDir:=eVectDir * LorMat * eVectDir;
  eVectBasDir:=eVectBas + eVectDir;
  eNormBasDir:=eVectBasDir * LorMat * eVectBasDir;
  GetMidVal:=function(TheLow, TheUpp)
      local eFrac, TargetLow, TargetUpp, TheSeq, eVal;
      eFrac:=(TheLow + TheUpp) / 2;
      TargetLow:=(2*TheLow + TheUpp) / 3;
      TargetUpp:=(TheLow + 2*TheUpp) / 3;
      TheSeq:=GetSequenceContinuousFractionApproximant(eFrac);
#      Print("TheSeq=", TheSeq, "\n");
      for eVal in TheSeq
      do
          if TargetLow <= eVal and eVal <= TargetUpp then
              return eVal;
          fi;
      od;
      Error("Should never reach that stage");
  end;
  while(true)
  do
#    Print("TheLowerBound=", TheLowerBound, " TheUpperBound=", TheUpperBound, "\n");
    TheMidVal:=GetMidVal(TheLowerBound, TheUpperBound);
    eNSPtest:=eNSPbas + TheMidVal * eNSPdir;
    eVectTest:=RemoveFraction(eNSPtest{[2..n+1]}*Inverse(LorMat));
    eNormTest:=eVectTest*LorMat*eVectTest;
    MaxScal:=CritSet[1]*LorMat*eVectTest;
#    Print("TheMidVal=", TheMidVal, " eNormTest=", eNormTest, " MaxScal=", MaxScal, " eNSPtest=", eNSPtest, " eVectTest=", eVectTest, "\n");
#    Print("TheMidVal=", TheMidVal, " eNormTest=", eNormTest, " MaxScal=", MaxScal, " eNSPtest=", eNSPtest, "\n");
#    Print("eNormBas=", eNormBas, " eNormDir=", eNormDir, " eNormBasDir=", eNormBasDir, "\n");
    if eNormTest <= 0 or MaxScal <= 0 then
      TheUpperBound:=TheMidVal;
    else
      ListTotal:=LORENTZ_FindPositiveVectors(LorMat, eVectTest, MaxScal, TheOption, OnlyShortest);
#      Print("ListTotal=", ListTotal, "\n");
#      Print("  |ListTotal|=", Length(ListTotal), "\n");
      if IsSubset(CritSet, Set(ListTotal)) and Length(CritSet) > Length(ListTotal) then
        Error("Bug: if included, it should be equal");
      fi;
      if Set(ListTotal)=Set(CritSet) then
        TheLowerBound:=TheMidVal;
      else
        if IsSubset(Set(ListTotal), CritSet) then
#          Print("EXIT 1 |ListTotal|=", Length(ListTotal), " NSPtest=", eNSPtest, " MaxScal=", MaxScal, "\n");
          return rec(ListTotal:=ListTotal, eNSPtest:=eNSPtest, eVectTest:=eVectTest, MaxScal:=MaxScal);
        else
          break;
        fi;
      fi;
    fi;
#    if n_iter > 10 then
#        Print("n_iter=", n_iter, " is too large\n");
#        Error("Stop here");
#    fi;
    n_iter:=n_iter+1;
  od;
#  Print("Going to the second scheme\n");
  while(true)
  do
    eVect:=ListTotal[1];
    eVert:=Concatenation([1], eVect);
    ListIsoTest:=Concatenation(CritSet, [eVect]);
    EXTtest:=List(ListIsoTest, x->Concatenation([1], x));
    aShift:=-( eNSPbas*eVert) / ( eNSPdir*eVert );
    eNSPtest:=eNSPbas + aShift*eNSPdir;
#    Print(" aShift=", aShift, "\n");
    eVectTest:=RemoveFraction(eNSPtest{[2..n+1]}*Inverse(LorMat));
    MaxScal:=CritSet[1]*LorMat*eVectTest;
    ListTotal:=LORENTZ_FindPositiveVectors(LorMat, eVectTest, MaxScal, TheOption, OnlyShortest);
    if IsSubset(Set(ListTotal), Set(ListIsoTest)) then
#      Print("EXIT 2 |ListTotal|=", Length(ListTotal), " NSPtest=", eNSPtest, " MaxScal=", MaxScal, "\n");
      return rec(ListTotal:=ListTotal, eNSPtest:=eNSPtest, eVectTest:=eVectTest, MaxScal:=MaxScal);
    fi;
  od;
end;


LORENTZ_GetShortPositiveVector:=function(LorMat)
    local L1_Norm, n, GraverBasis, k, LVect, eVect, fVect, eSet, DirectImprovement, ePerturb, n_no_improv, CurrNorm, CurrVect, LorMat_Pert, TheVect, eNorm1, eNorm, norm1, norm2, uVect, ListCurrVect, FuncInsert, i, GetAlpha;
    n:=Length(LorMat);
    L1_Norm:=function(eMat)
        return Sum(List(eMat, x->Sum(List(x, AbsInt))));
    end;
    GraverBasis:=[];
    for k in [1..2]
    do
        LVect:=BuildSet(k, [-1,1]);
        for eSet in Combinations([1..n], k)
        do
            eVect:=ListWithIdenticalEntries(n,0);
            for fVect in LVect
            do
                for i in [1..k]
                do
                    eVect[eSet[i]] := fVect[i];
                od;
                Add(GraverBasis, ShallowCopy(eVect));
            od;
        od;
    od;
    GetAlpha:=function(TheVect, DirVect)
        local Norm1, Norm2, NewVect, alpha;
        Norm1:=TheVect * LorMat * TheVect;
        alpha:=0;
        while(true)
        do
            NewVect:=TheVect + (alpha+1) * DirVect;
            Norm2:=NewVect * LorMat * NewVect;
#            Print("Norm1=", Norm1, " Norm2=", Norm2, " alpha=", alpha, "\n");
            if Norm2 < Norm1 and Norm2 > 0 then
                Norm1:=Norm2;
                alpha:=alpha+1;
            else
                return alpha;
            fi;
        od;
    end;
    DirectImprovement:=function(TheVect)
        local ImpNorm, ImpVect, NoImprov, eVectBasis, NewVect, NewNorm, n_chg, alpha;
        ImpVect:=TheVect;
        ImpNorm:=TheVect * LorMat * TheVect;
        while(true)
        do
            n_chg:=0;
            for eVectBasis in GraverBasis
            do
                alpha:=GetAlpha(ImpVect, eVectBasis);
                ImpVect:=ImpVect + alpha * eVectBasis;
#                Print("eVectBasis=", eVectBasis, " alpha=", alpha, "\n");
                n_chg:=n_chg + alpha;
            od;
#            Print("n_chg=", n_chg, "\n");
            if n_chg = 0 then
                return ImpVect;
            fi;
        od;
    end;
    ePerturb:=IdentityMat(n);
    n_no_improv:=0;
    CurrNorm:=infinity;
    CurrVect:="unset";
    ListCurrVect:=[];
    FuncInsert:=function(TheVect)
        local ListNorm, TheMin, ListIdx;
        Add(ListCurrVect, TheVect);
        ListNorm:=List(ListCurrVect, x->x * LorMat * x);
        TheMin:=Minimum(ListNorm);
        ListIdx:=Filtered([1..Length(ListNorm)], x->ListNorm[x] <= 3*TheMin);
        ListCurrVect:=ListCurrVect{ListIdx};
    end;
    while(true)
    do
        LorMat_Pert:=ePerturb * LorMat * TransposedMat(ePerturb);
        uVect:=FindNegativeVect(-LorMat_Pert);
        eNorm1:=uVect * LorMat_Pert * uVect;
        TheVect:=DirectImprovement(uVect * ePerturb);
        eNorm:=TheVect * LorMat * TheVect;
        if eNorm <= 0 then
            Print("n_no_improv=", n_no_improv, "  eNorm=", eNorm, " TheVect=", TheVect, " CurrNorm=", CurrNorm, "\n");
            Error("eNorm should be strictly positive");
        fi;
        if eNorm < CurrNorm then
            n_no_improv:=0;
            CurrVect:=TheVect;
            CurrNorm:=eNorm;
        else
            n_no_improv:=n_no_improv+1;
        fi;
        if n_no_improv=500 or CurrNorm < 10 then
            if CurrVect="unset" then
                Error("CurrVect is unset");
            fi;
            return CurrVect;
        fi;
        norm1:=L1_Norm(LorMat_Pert);
        norm2:=L1_Norm(LorMat);
        if norm1 > 10000 * norm2 then
            ePerturb:=IdentityMat(n);
        fi;
        ePerturb:=ePerturb * GetRandomMatrixPerturbation(n);
    od;
end;


# Given a Lorentzian form, find:
# ListTotal: vectors that correspond to a perfect form
# eNSPtest: A vector v such that v*w = 0 for all w in ListTotal
# CentralVect: A vectors v of positive norm that correspond essentially
#    to one of the selected cones.
LORENTZ_GetOnePerfect:=function(LorMat, TheOption)
    local eRec, n, pos, eVect, eScal, CritSet, eNSPbas, EXT, NSP, eVEctB, eNSPdir, eRecB, rnk, eVectB, CentralVect, IsInCone, ListDir, GetOneOutsideRay, Viso, CheckVectorEXT;
    n:=Length(LorMat);
#    Print("Beginning of LORENTZ_GetOnePerfect, TheOption=", TheOption, "LorMat=\n");
#    PrintArray(LorMat);
    if LORENTZ_IsLorentzian(LorMat)=false then
        Error("LorMat should be Lorentzian");
    fi;
    GetOneOutsideRay:=function(SpannBasis, TheSet)
        local TheMat, n_dim, uVect, RetVect, ListScal, eScal, TheNorm, ePerturb, TheMatPerturb;
        TheMat:=SpannBasis * LorMat * TransposedMat(SpannBasis);
#        Print("From SpannBasis, TheMat=\n");
#        PrintArray(TheMat);
        if DiagonalizeSymmetricMatrix(TheMat).nbMinus=0 then
            Error("We should have a negative entry in the matrix");
        fi;
        n_dim:=Length(TheMat);
        ePerturb:=IdentityMat(n_dim);
        while(true)
        do
            TheMatPerturb:=ePerturb * TheMat * TransposedMat(ePerturb);
            uVect:=FindNegativeVect(TheMatPerturb);
            RetVect:=uVect * ePerturb * SpannBasis;
            TheNorm:=RetVect * LorMat * RetVect;
#            Print("TheNorm=", TheNorm, "\n");
            if TheNorm >= 0 then
                Error("The vector should be outside of the cone and so have negative norm");
            fi;
            ListScal:=List(TheSet, x->x*LorMat*RetVect);
            if Length(Set(ListScal))<>1 then
                Error("The scalar products should all be the same");
            fi;
            eScal:=ListScal[1];
            eNSPdir:=Concatenation([-eScal], LorMat*RetVect);
            # We need to check for the bad scenario of obtaining an isotropic space 
            if GetUpperBound(LorMat, eNSPbas, eNSPdir)<>fail then
                return eNSPdir;
            fi;
            ePerturb:=ePerturb * GetRandomMatrixPerturbation(n_dim);
        od;
    end;
    CheckVectorEXT:=function(EXT, eVect)
        if First(EXT, x->x*eVect<>0)<>fail then
            Error("Please debug from here Lorentzian CheckVectorEXT");
        fi;
    end;
    CentralVect:=LORENTZ_GetShortPositiveVector(LorMat);
#    Print("CentralVect=", CentralVect, "\n");
    if n > 4 and false then # In that case an isotropic vector
#        Print("Running isotropic code\n");
        Viso:=INDEF_FindIsotropic(LorMat);
        eScal:=Viso * LorMat * CentralVect;
#        Print("eScal=", eScal, "\n");
        if eScal < 0 then
            Viso:=-Viso;
        fi;
        CritSet:=[Viso];
    else
#        Print("Running classic SearchInitialVector code\n");
        CritSet:=LORENTZ_SearchInitialVector(LorMat, CentralVect, TheOption);
    fi;
    eScal:=CritSet[1]*LorMat*CentralVect;
#    Print("|CritSet|=", Length(CritSet), " eScal=", eScal, " l_norm=", List(CritSet, x->x*LorMat*x), "\n");
    eNSPbas:=Concatenation([-eScal], CentralVect*LorMat);
    while(true)
    do
        rnk:=RankMat(CritSet);
#        Print("LORENTZ_GetOnePerfect: rnk=", rnk, " ListTotal=", Set(CritSet), "\n");
        if rnk = n then
            LORENTZ_CheckCorrectnessVectorFamily(CritSet);
            return rec(ListTotal:=CritSet, eNSPtest:=eNSPbas, eVectTest:=CentralVect);
        fi;
        EXT:=List(CritSet, x->Concatenation([1], x));
        CheckVectorEXT(EXT, eNSPbas);
        NSP:=NullspaceMat(TransposedMat(EXT));
        if Length(NSP)=0 then
            Error("NSP should be non-empty");
        fi;
#        Print("LORENTZ_GetOnePerfect: |NSP|=", Length(NSP), "\n");
        ListDir:=List(NSP, x->RemoveFraction(x{[2..n+1]}*Inverse(LorMat)));
        eNSPdir:=GetOneOutsideRay(ListDir, CritSet);
        CheckVectorEXT(EXT, eNSPdir);
#        Print("LORENTZ_GetOnePerfect: eNSPdir=", eNSPdir, "\n");
        eRecB:=LORENTZ_Kernel_Flipping(LorMat, CritSet, eNSPbas, eNSPdir, TheOption);
        CritSet:=eRecB.ListTotal;
        eNSPbas:=RemoveFraction(eRecB.eNSPtest);
    od;
end;


# The function computes the flipping.
# It is essentially just a syntactic sugar of the LORENTZ_Kernel_Flipping
LORENTZ_DoFlipping:=function(LorMat, ListIso, eInc, TheOption)
  local n, EXT, NSP, eVert, eEXT, eNSP, NSPb, eNSPb, iShift, CritSet, eNSPtest, eVectTest, MaxScal, ListTotal, EXTtest, ListIsoTest, NSPtest, eVectB, eNormTest, eNSPdir, eNSPbas, TheFlip;
  if LORENTZ_IsLorentzian(LorMat)=false then
    Error("LorMat should be Lorentzian");
  fi;
  n:=Length(LorMat);
  EXT:=List(ListIso, x->Concatenation([1], x));
  eVert:=Difference([1..Length(EXT)], eInc)[1];
  eEXT:=EXT[eVert];
  NSP:=NullspaceMat(TransposedMat(ListIso{eInc}));
  if Length(NSP)<>1 then
    Error("The array NSP should have size 1");
  fi;
  if NSP[1]*ListIso[eVert]<0 then
    eNSPdir:=RemoveFraction(Concatenation([0], -NSP[1]));
  else
    eNSPdir:=RemoveFraction(Concatenation([0],  NSP[1]));
  fi;
  NSPb:=NullspaceMat(TransposedMat(EXT));
  eVectB:=NSPb[1]{[2..n+1]};
  if eVectB*ListIso[eVert]>0 then
    eNSPbas:=RemoveFraction(NSPb[1]);
  else
    eNSPbas:=RemoveFraction(-NSPb[1]);
  fi;
  CritSet:=Set(ListIso{eInc});
  TheFlip:=LORENTZ_Kernel_Flipping(LorMat, CritSet, eNSPbas, eNSPdir, TheOption).ListTotal;
  LORENTZ_CheckCorrectnessVectorFamily(TheFlip);
  return TheFlip;
end;

LORENTZ_Invariant:=function(LorMat, eFamEXT)
  return GetScalarMatrixInvariant_Polytope_AddMat(eFamEXT, [LorMat]);
end;


LORENTZ_EnumeratePerfect_DelaunayScheme:=function(LorMat, RecInput)
  local n, FuncStabilizer, FuncIsomorphy, FuncInvariant, WorkingDim, IsBankSave, IsRespawn, BF, FindDelaunayPolytope, FindAdjacentDelaunay, KillingDelaunay, KillingAdjacency, DataLattice, DataPolyhedral, TheReply, DelaunayDatabase, EXT, eStab, FuncStabilizerDelaunay, FuncIsomorphismDelaunay, MaximalMemorySave, ListFamily, iOrb, TheRec, TheOption, ListVertAnsatz, ListAnsatzInfo, TestNeedMoreSymmetry, testSym, ListAdj, eInv;
  n:=Length(LorMat)-1;
  FuncStabilizer:=LinPolytope_Automorphism;
  FuncIsomorphy:=LinPolytope_Isomorphism;
  FuncInvariant:=LinPolytope_Invariant;
  TheOption:="isotrop";
  if IsBound(RecInput.TheOption) then
    TheOption:=RecInput.TheOption;
  fi;
  KillingDelaunay:=function(EXT, eInv)
    return false;
  end;
  if IsBound(RecInput.KillingDelaunay) then
    KillingDelaunay:=RecInput.KillingDelaunay;
  fi;
  WorkingDim:=Length(LorMat);
  IsBankSave:=function(OrdGRP, EXT, TheDepth)
    if TheDepth=0 then
      return false;
    fi;
    if Length(EXT)<=WorkingDim+5 then
      return false;
    fi;
    return true;
  end;
  if IsBound(RecInput.IsBankSave) then
    IsBankSave:=RecInput.IsBankSave;
  fi;
  IsRespawn:=function(OrdGRP, EXT, TheDepth)
    if OrdGRP>=50000 and TheDepth<=2 then
      return true;
    fi;
    if OrdGRP<100 then
      return false;
    fi;
    if Length(EXT)<WorkingDim+20 then
      return false;
    fi;
    if TheDepth=2 then
      return false;
    fi;
    return true;
  end;
  if IsBound(RecInput.IsRespawn) then
    IsRespawn:=RecInput.IsRespawn;
  fi;
  BF:=BankRecording(FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism(500));
  TestNeedMoreSymmetry:=function(EXT)
    if Length(EXT) > RankMat(EXT) + 4 then
      return true;
    else
      return false;
    fi;
  end;
  if IsBound(RecInput.TestNeedMoreSymmetry) then
      testSym:=RecInput.TestNeedMoreSymmetry;
  else
      testSym:=TestNeedMoreSymmetry;
  fi;
  DataPolyhedral:=rec(IsBankSave:=IsBankSave,
        TheDepth:=0,
        IsRespawn:=IsRespawn,
        GetInitialRays:=GetInitialRays_LinProg,
        TestNeedMoreSymmetry:=testSym,
        FuncStabilizer:=FuncStabilizer,
        FuncIsomorphy:=FuncIsomorphy,
        FuncInvariant:=FuncInvariant,
        DualDescriptionFunction:=__DualDescriptionLRS_Reduction,
        GroupFormalism:=OnSetsGroupFormalism(500));
  #
  # The geometrical part
  #
  KillingAdjacency:=function(EXT1, EXT2)
    return false;
  end;
  ListAnsatzInfo:=[];
  if IsBound(RecInput.ListAnsatzInfo) then
    ListAnsatzInfo:=RecInput.ListAnsatzInfo;
  fi;
  ListVertAnsatz:=List(ListAnsatzInfo, x->x.nbVert);
  FuncStabilizerDelaunay:=function(eRec, EXT)
    local pos, RecRet;
    RecRet:=LORENTZ_ComputeStabilizer(LorMat, EXT);
    RecRet.PermutationStabilizer:=RecRet.GRP_int;
    return RecRet;
  end;
  FuncIsomorphismDelaunay:=function(eRec, EXT1, EXT2, eStab1)
    local pos;
    if Length(EXT1)<>Length(EXT2) then
      return false;
    fi;
    return LORENTZ_TestEquivalence(LorMat, EXT1, EXT2);
  end;
  FindDelaunayPolytope:=function()
    return LORENTZ_GetOnePerfect(LorMat, TheOption).ListTotal;
  end;
  FuncInvariant:=function(eRec, EXT)
    return LORENTZ_Invariant(LorMat, EXT);
  end;
  FindAdjacentDelaunay:=function(EXT, eOrb)
    local ListVect;
    ListVect:=List(EXT, x->x{[2..n+1]});
    return LORENTZ_DoFlipping(LorMat, EXT, eOrb, TheOption);
  end;
  DataLattice:=rec(n:=n,
                   KillingDelaunay:=KillingDelaunay,
                   KillingAdjacency:=KillingAdjacency,
                   FindDelaunayPolytope:=FindDelaunayPolytope,
                   FindAdjacentDelaunay:=FindAdjacentDelaunay,
                   FuncInvariant:=FuncInvariant,
		   FuncIsomorphismDelaunay:=FuncIsomorphismDelaunay,
                   FuncStabilizerDelaunay:=FuncStabilizerDelaunay);
  #
  # The saving business part
  #
  DelaunayDatabase:=DelaunayDatabaseManagement();
  TheReply:=ComputeDelaunayDecomposition(DataLattice, DataPolyhedral, DelaunayDatabase);
  if TheReply<>"all was ok" then
    return TheReply;
  fi;
  ListFamily:=[];
  for iOrb in [1..DelaunayDatabase.FuncDelaunayGetNumber()]
  do
    EXT:=DelaunayDatabase.FuncDelaunayGetEXT(iOrb);
    eInv:=DelaunayDatabase.FuncDelaunayGetINV(iOrb);
    eStab:=DelaunayDatabase.FuncDelaunayGetGroup(iOrb);
    ListAdj:=DelaunayDatabase.FuncDelaunayGetAdjacencies(iOrb);
    TheRec:=rec(EXT:=EXT, eInv:=eInv, eStab:=eStab, ListAdj:=ListAdj);
    Add(ListFamily, TheRec);
  od;
  return ListFamily;
end;


LORENTZ_InvariantOfLorentzianForm:=function(LorMat)
  return DeterminantMat(LorMat);
end;


#
# It is very heavy handed (compute all perfect forms ...)
# but it works
LORENTZ_TestIsomorphismLorentzianMatrices:=function(LorMat1, LorMat2)
  local eInv1, eInv2, TheOption, RecInput, eFamEXT2, eInvPerfect, KillingDelaunay, TheResult;
  eInv1:=LORENTZ_InvariantOfLorentzianForm(LorMat1);
  eInv2:=LORENTZ_InvariantOfLorentzianForm(LorMat2);
  if eInv1<>eInv2 then
    return false;
  fi;
  TheOption:="total";
  eFamEXT2:=LORENTZ_GetOnePerfect(LorMat2, TheOption).ListTotal;
  eInvPerfect:=LORENTZ_Invariant(LorMat2, eFamEXT2);
#  Print("|eFamEXT2|=", Length(eFamEXT2), " eInvPerfect=", eInvPerfect, "\n");
  #
  KillingDelaunay:=function(eFamEXT1, eInv)
      local test, testInv;
#      Print("|eFamEXT1|=", Length(eFamEXT1), " eInv=", eInv, "\n");
      if eInv=eInvPerfect then
          test:=LORENTZ_TestEquivalence_General(LorMat1, LorMat2, eFamEXT1, eFamEXT2);
          if test<>false then
              testInv:=Inverse(test);
              if testInv*LorMat1*TransposedMat(testInv)<>LorMat2 then
                  Error("The program still has some bugs");
              fi;
              return rec(TheEquiv:=testInv); # Need to encapsulate the matrix into a record, otherwise it matches the IsList
          fi;
      fi;
      return false;
  end;
  RecInput:=rec(TheOption:=TheOption, KillingDelaunay:=KillingDelaunay);
  TheResult:=LORENTZ_EnumeratePerfect_DelaunayScheme(LorMat1, RecInput);
  if IsList(TheResult) then # TheResult is then the list of perfect forms
      return false; # None of the perfect forms matched
  fi;
  return TheResult.TheEquiv; # Has to be a record
end;


LORENTZ_OrthogonalGroup_Kernel:=function(LorMat, ListFamily)
    local n, ListGenerators, ePerf, eGen, eAdj, SetGen;
    n:=Length(LorMat);
    ListGenerators:=[ - IdentityMat(n)]; # We need a transformation that flips the cone. The other below dont.
    for ePerf in ListFamily
    do
        for eGen in GeneratorsOfGroup(ePerf.eStab.GRPintMatr)
        do
            Add(ListGenerators, eGen);
        od;
        for eAdj in ePerf.ListAdj
        do
            Add(ListGenerators, eAdj.eBigMat);
        od;
    od;
    SetGen:=Set(ListGenerators);
    for eGen in SetGen
    do
        if eGen * LorMat * TransposedMat(eGen) <> LorMat then
            Error("The matrix does not preserve the quadratic form");
        fi;
        if IsIntegralMat(eGen)=false then
            Error("The matrix is not integral");
        fi;
    od;
    return Group(SetGen);
end;


LORENTZ_OrthogonalGroup:=function(LorMat)
    local RecInput, ListFamily;
    RecInput:=rec(TheOption:="total");
    ListFamily:=LORENTZ_EnumeratePerfect_DelaunayScheme(LorMat, RecInput);
    return LORENTZ_OrthogonalGroup_Kernel(LorMat, ListFamily);
end;


LORENTZ_CandidateIsotropicVectors_Kernel:=function(LorMat, ListFamily)
    local n, ListCand, ePerf, GRPperm, EXT, O, eO, eVectCand, eNorm;
    n:=Length(LorMat);
    ListCand:=[];
    for ePerf in ListFamily
    do
        GRPperm:=ePerf.eStab.PermutationStabilizer;
        EXT:=ePerf.EXT;
        O:=Orbits(GRPperm, [1..Length(EXT)], OnPoints);
        for eO in O
        do
            eVectCand:=EXT[eO[1]];
            eNorm:=eVectCand * LorMat * eVectCand;
            if eNorm=0 then
                Add(ListCand, eVectCand);
            fi;
        od;
    od;
    return ListCand;
end;
