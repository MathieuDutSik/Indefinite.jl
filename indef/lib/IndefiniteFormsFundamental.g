IndefiniteReduction:=function(M)
    local M_oscar, TheResult_oscar, Mred, B;
    M_oscar:=MatrixToOscar(M);
    TheResult_oscar:=Oscar.IndefiniteReduction(M_oscar);
    Mred:=ReadOscarMatrix(TheResult_oscar[1]);
    B:=ReadOscarMatrix(TheResult_oscar[2]);
    return rec(Mred:=Mred, B:=B);
end;


GetRandomMatrixPerturbation:=function(n)
    local choice, ePerm, eMat, i, j;
    choice:=Random([1..2]);
    if choice = 1 then
        ePerm:=Random(SymmetricGroup(n));
        eMat:=NullMat(n,n);
        for i in [1..n]
        do
            eMat[i][OnPoints(i, ePerm)]:=Random([-1,1]);
        od;
        return eMat;
    fi;
    if choice = 2 then
        eMat:=IdentityMat(n);
        i:=Random([1..n]);
        j:=Random([1..n]);
        eMat[i][j]:=Random([-1,1]);
        return eMat;
    fi;
end;


INDEF_GetIndefinitesubset:=function(M)
    local n, eSet, Mred, RecDiag;
    n:=Length(M);
    if n <= 9 then
        return rec(subset:=[1..n], Mred:=M, look_kernel:=false);
    fi;
    for eSet in Combinations([1..n],9)
    do
        Mred:=List(M{eSet}, x->x{eSet});
        RecDiag:=DiagonalizeSymmetricMatrix(Mred);
        if RecDiag.nbZero>0 or (RecDiag.nbMinus > 0 and RecDiag.nbPlus > 0) then
            return rec(subset:=eSet, Mred:=Mred, look_kernel:=RecDiag.nbZero > 0);
        fi;
    od;
    Print("Failed to find a subset");
end;


INDEF_FindIsotropic:=function(M)
    local dim, M_oscar, q, reply, v, v_list, eList, eVect, retVect;
    dim:=Length(M);
    Print("M=", M, "\n");
    M_oscar:=MatrixToOscar(M);
    Print("M_oscar=", M_oscar, "\n");
    q := Oscar.quadratic_space(Oscar.QQ, M_oscar);
    reply := Oscar.is_isotropic_with_vector(q);
    Print("reply=", reply, "\n");
    v:=reply[2];
    Print("v=", v, "\n");
    v_list:=JuliaToGAP(IsList, v);
    eVect:=List(v_list, Oscar.GAP.julia_to_gap);
    Print("eVect=", eVect, "\n");
    retVect:=RemoveFraction(eVect);
    return retVect;
end;


INDEF_FORM_IsEven:=function(Qmat)
    local n;
    if IsIntegralMat(Qmat)=false then
        return false;
    fi;
    n:=Length(Qmat);
    return First([1..n], x->(Qmat[x][x] mod 2)=1)=fail;
end;


INDEF_FORM_Invariant:=function(Qmat)
    local n, NSP, TheCompl, GramRed, eDet, DiagInfo, nbPlus, nbMinus, nbZero, IsEven;
    if IsIntegralMat(Qmat)=false then
        Error("We consider only integral matrices");
    fi;
    n:=Length(Qmat);
    NSP:=NullspaceIntMat(Qmat);
    TheCompl:=SubspaceCompletionInt(NSP, n);
    GramRed:=TheCompl * Qmat * TransposedMat(TheCompl);
    eDet:=DeterminantMat(GramRed);
    DiagInfo:=DiagonalizeSymmetricMatrix(GramRed);
    nbPlus:=DiagInfo.nbPlus;
    nbMinus:=DiagInfo.nbMinus;
    nbZero:=Length(NSP);
    IsEven:=INDEF_FORM_IsEven(Qmat);
    return rec(n:=n, eDet:=eDet, nbPlus:=nbPlus, nbMinus:=nbMinus, nbZero:=nbZero, IsEven:=IsEven);
end;


INDEF_FORM_InvariantVector:=function(Qmat, v)
    local eNorm, eProd, divisor, index, NSP, GramRed, typeInv, TheRank;
    TheRank:=RankMat(Qmat);
    eNorm:=v * Qmat * v;
    eProd:=v * Qmat;
    divisor:=1 / RemoveFractionPlusCoef(v).TheMult;
    index:=GcdVector(eProd).TheGcd;
    NSP:=NullspaceIntMat(TransposedMat([v * Qmat]));
    GramRed:=NSP * Qmat * TransposedMat(NSP);
    typeInv:=INDEF_FORM_Invariant(GramRed);
    return rec(eRank:=TheRank,
               eNorm:=eNorm,
               index:=index,
               divisor:=divisor,
               typeInv:=typeInv);
end;


INDEF_FORM_Invariant_IsotropicKplane_Raw:=function(Qmat, ePlane)
    local k, NSP, dimNSP, ePlaneB, eV, eSol, ComplBasisInNSP, NSP_sub, QmatRed, eInv1, eInv2;
    k:=Length(ePlane);
    NSP:=NullspaceIntMat(TransposedMat(ePlane * Qmat));
    dimNSP:=Length(NSP);
    ePlaneB:=[];
    for eV in ePlane
    do
        eSol:=SolutionIntMat(NSP, eV);
        if eSol=fail then
            Error("eV should belong to the space by the virtue of being isotropic");
        fi;
        Add(ePlaneB, eSol);
    od;
    ComplBasisInNSP:=SubspaceCompletionInt(ePlaneB, dimNSP);
    NSP_sub:=ComplBasisInNSP * NSP;
    QmatRed:=NSP_sub * Qmat * TransposedMat(NSP_sub);
    eInv1:=INDEF_FORM_Invariant(Qmat);
    if DeterminantMat(QmatRed)=0 then
        Error("QmatRed should be non-degenerate");
    fi;
    eInv2:=INDEF_FORM_Invariant(QmatRed);
    return rec(k:=k, eInv1:=eInv1, eInv2:=eInv2);
end;


TestEqualitySpace:=function(Space1, Space2)
    local eV1, eV2;
    if Length(Space1)<>Length(Space2) then
        return false;
    fi;
    for eV1 in Space1
    do
        if SolutionIntMat(Space2, eV1)=fail then
            return false;
        fi;
    od;
    for eV2 in Space2
    do
        if SolutionIntMat(Space1, eV2)=fail then
            return false;
        fi;
    od;
    return true;
end;
