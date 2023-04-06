FileConvertPariIsotropOutput:=Filename(DirectoriesPackagePrograms("indefinite"),"ConvertPariIsotrop");
FileIndefiniteReduction:=Filename(DirectoriesPackagePrograms("indefinite"),"IndefiniteReduction");


IndefiniteReduction:=function(M)
    local FileI, FileO, FileE, n, output, i, j, eVal, eCommand, TheReply;
    FileI:=Filename(POLYHEDRAL_tmpdir,"Indefinite.input");
    FileO:=Filename(POLYHEDRAL_tmpdir,"Indefinite.output");
    FileE:=Filename(POLYHEDRAL_tmpdir,"Indefinite.error");
#    Print("FileI=", FileI, "\n");
#    Print("FileO=", FileO, "\n");
#    Print("FileE=", FileE, "\n");
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    RemoveFileIfExist(FileE);
    n:=Length(M);
    output:=OutputTextFile(FileI, true);
    AppendTo(output, n, " ", n, "\n");
    for i in [1..n]
    do
        for j in [1..n]
        do
            eVal:=M[i][j];
            AppendTo(output, " ", eVal);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
    #
    eCommand:=Concatenation(FileIndefiniteReduction, " ", FileI, " ", FileO, " 2> ", FileE);
#    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    TheReply:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    RemoveFileIfExist(FileE);
    return TheReply;
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


INDEF_FindIsotropic_Kernel:=function(M)
    local n, FileInput, FileOutput, FileRead, FileErr, RecGRPcone, output, TheRec, TheCommand1, TheCommand2, eVect, ListCorrect, ListNorm, MinNorm, pos, eRes, FullVect, ClearFiles;
    n:=Length(M);
    FileInput :=Filename(POLYHEDRAL_tmpdir,"Isotrop.pari");
    FileOutput:=Filename(POLYHEDRAL_tmpdir,"Isotrop.out");
    FileRead:=Filename(POLYHEDRAL_tmpdir,"Isotrop.read");
    FileErr:=Filename(POLYHEDRAL_tmpdir,"Isotrop.err");
    RemoveFileIfExist(FileInput);
    #
    TheRec:=INDEF_GetIndefinitesubset(M);
    #
    output:=OutputTextFile(FileInput, true);
    AppendTo(output, "M=");
    PARI_PrintMatrix(output, TheRec.Mred);
    AppendTo(output, "\n");
    AppendTo(output, "eR = iferr(qfsolve(M),E,[],errname(E)==\"e_IMPL\")\n");
    AppendTo(output, "print(eR)\n");
    AppendTo(output, "quit\n");
    CloseStream(output);
    #
    TheCommand1:=Concatenation("gp ", FileInput, " > ", FileOutput, " 2> ", FileErr);
#    Print("TheCommand1=", TheCommand1, "\n");
    Exec(TheCommand1);
#    Print("TheCommand1 has been executed\n");
    #
    TheCommand2:=Concatenation(FileConvertPariIsotropOutput, " ", FileOutput, " ", FileRead);
#    Print("TheCommand2=", TheCommand2, "\n");
    Exec(TheCommand2);
    #
    eRes:=ReadAsFunction(FileRead)();
    ClearFiles:=function()
        RemoveFileIfExist(FileInput);
        RemoveFileIfExist(FileOutput);
        RemoveFileIfExist(FileRead);
        RemoveFileIfExist(FileErr);
    end;
    ListCorrect:=[];
    for eVect in eRes
    do
        if Length(eVect) = Length(TheRec.Mred) then
            FullVect:=ListWithIdenticalEntries(Length(M),0);
            FullVect{TheRec.subset}:=eVect;
            if FullVect * M * FullVect = 0 then
                Add(ListCorrect, FullVect);
            fi;
        fi;
    od;
    if Length(ListCorrect) = 0 then
        Print("gp seems to have failed. Vector is not isotrop");
#        SaveDebugInfo("FailureIsotropSearch", M);
        ClearFiles();
        return fail;
    fi;
    ListNorm:=List(ListCorrect, Norm_L1);
    MinNorm:=Minimum(ListNorm);
    pos:=Position(ListNorm, MinNorm);
    ClearFiles();
    return ListCorrect[pos];
end;


INDEF_FindIsotropic:=function(M)
    local n_iter, n, ThePerturb, ePerturb, Mtest, TheReply, OneIso, eScal;
    n_iter:=0;
    n:=Length(M);
    ThePerturb:=IdentityMat(n);
    while(true)
    do
        ePerturb:=GetRandomMatrixPerturbation(n);
        ThePerturb:=ePerturb * ThePerturb;
        Mtest:=ThePerturb * M * TransposedMat(ThePerturb);
        TheReply:=INDEF_FindIsotropic_Kernel(Mtest);
        if TheReply<>fail then
            OneIso:=TheReply * ThePerturb;
            eScal:=OneIso * M * OneIso;
            if eScal<>0 then
                Error("The eScal should be zero");
            fi;
            return OneIso;
        fi;
        n_iter:=n_iter+1;
        Print("Retrying at n_iter=", n_iter, "\n");
    od;
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
