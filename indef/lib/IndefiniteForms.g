

INDEF_FORM_GetAttackScheme:=function(Qmat)
    local DiagInfo, nbPlus, nbMinus;
    DiagInfo:=DiagonalizeSymmetricMatrix(Qmat);
    nbPlus:=DiagInfo.nbPlus;
    nbMinus:=DiagInfo.nbMinus;
    Print("nbPlus=", nbPlus, " nbMinus=", nbMinus, "\n");
    if nbMinus <= nbPlus then
        return rec(h:=nbMinus, mat:=-Qmat);
    fi;
    return rec(h:=nbPlus, mat:=Qmat);
end;


# The Eichler transvection are defined for an even lattice (so integral and even norms of vectors.
# For such a lattice L* is defined as Q^{-1} Z^n
# The formula describing them is
# E_{f,x}(y) = y + (y,x) f - (x,x) / 2 (y,f) f - (y,f) x
# The Eichler transvection preserve the lattice L because (x,x) / 2 in Z, (y,x) in Z and (y,f) in Z.
# Also if y in L* since the above is still true, we have that E_{f,x} preserve L* and the image of
# the transformation in L*/L is the identity.
#
# Test of composition property:
# We have
# E_{f,y}(z) = z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y
# E_{f,x}E_{f,y}(z) =
#   (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y)
#   + (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y, x) f
#   - (x,x) / 2 (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y, f) f
#   - (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y, f) x
# = z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y
#   + (z - (z,f) y, x) f
#   - (x,x) / 2 (z,f) f
#   - (z,f) x
# = z + (z,x + y) f - (z,f) (x+y) - (y,y)/2 (z,f) f - (x,x)/2 (z,f) f
#     -(z,f) (y,x) f
# = z + (z,x + y) f - (z,f) (x+y) - (x+y,x+y)/2 (z,f) f
INDEF_FORM_Eichler_Transvection:=function(Qmat, f, x)
    local n, fNorm, scal, xNorm, RetMat, y, eImg, scal_yx, scal_yf;
    if INDEF_FORM_IsEven(Qmat)=false then
        Error("The lattice Qmat should be even in order to define the Eichler transvection");
    fi;
    n:=Length(Qmat);
    fNorm:=f * Qmat * f;
    scal:=f * Qmat * x;
    if fNorm<>0 or scal<>0 then
        Error("eNorm or scal are inconsistent");
    fi;
    xNorm:=x * Qmat * x;
    RetMat:=[];
    for y in IdentityMat(n)
    do
        eImg:=ListWithIdenticalEntries(n,0);
        eImg:=eImg + y;
        #
        scal_yx:=y * Qmat * x;
        eImg:=eImg + scal_yx * f;
        #
        scal_yf:=y * Qmat * f;
        eImg:=eImg - (xNorm/2) * scal_yf * f;
        #
        eImg:=eImg - scal_yf * x;
        #
        Add(RetMat, eImg);
    od;
    if RetMat * Qmat * TransposedMat(RetMat)<>Qmat then
        Error("The Eichler transvection is not an isometry");
    fi;
    if IsIntegralMat(RetMat)=false then
        Error("The Eichler transvection matrix should be integral");
    fi;
    return RetMat;
end;


SubBlock:=function(eMat, ListRow, ListCol)
    return List(eMat{ListRow}, x->x{ListCol});
end;


GetSquareDivisors:=function(X)
    local LFact, LColl, ListPrime, ListExpo, len, ListListProd, i, eVal, ListProd, u, ListSquareDiv, eEnt, eProd;
    LFact:=Factors(X);
    LColl:=Collected(LFact);
    ListPrime:=List(LColl, x->x[1]);
    ListExpo:=List(LColl, x->1 + LowerInteger(x[2]/2));
    len:=Length(ListExpo);
    ListListProd:=[];
    for i in [1..len]
    do
        eVal:=1;
        ListProd:=[];
        for u in [1..ListExpo[i]]
        do
            Add(ListProd, eVal);
            eVal:=eVal * ListPrime[i];
        od;
        Add(ListListProd, ListProd);
    od;
    ListSquareDiv:=[];
    for eEnt in BuildSetMultiple(ListListProd)
    do
        eProd:=Product(eEnt);
        Add(ListSquareDiv, eProd);
    od;
    return ListSquareDiv;
end;


INDEF_FORM_TestEquivalence_PosNeg:=function(Qmat1, Qmat2)
    local eRec1, eRec2, GetEquiv;
    eRec1:=DiagonalizeSymmetricMatrix(Qmat1);
    eRec2:=DiagonalizeSymmetricMatrix(Qmat1);
    if eRec1.strSign <> eRec2.strSign then
        return fail;
    fi;
    GetEquiv:=function(M1, M2)
        local test;
        test:=ArithmeticIsomorphism([M1], [M2]);
        if test=false then
            return fail;
        fi;
        return test;
    end;
    if eRec1.nbPlus > 0 then
        return GetEquiv(Qmat1, Qmat2);
    fi;
    return GetEquiv(-Qmat1, -Qmat2);
end;


INDEF_FORM_AutomorphismGroup_PosNeg:=function(Qmat)
    local eRec, GetEquiv;
    eRec:=DiagonalizeSymmetricMatrix(Qmat);
    if eRec.nbPlus > 0 then
        return ArithmeticAutomorphismGroup([Qmat]);
    fi;
    return ArithmeticAutomorphismGroup([-Qmat]);
end;


INDEF_FORM_GetOrbitRepresentative_PosNeg:=function(Qmat, X)
    local eRec, GetRepresentatives;
    eRec:=DiagonalizeSymmetricMatrix(Qmat);
    GetRepresentatives:=function(M, eNorm)
        local GRP, SHV1, SHV2, O, ListOrbit;
        if eNorm = 0 then
            return [];
        fi;
        GRP:=ArithmeticAutomorphismGroup([M]);
        SHV1:=ShortVectorDutourVersion(M, eNorm);
        SHV2:=Filtered(SHV1, x->x*M*x=eNorm);
        O:=Orbits(GRP, SHV2, OnPoints);
        ListOrbit:=List(O, x->x[1]);
        return ListOrbit;
    end;
    if eRec.nbPlus > 0 then
        return GetRepresentatives(Qmat, X);
    fi;
    return GetRepresentatives(Qmat, -X);
end;


# Based on paragraph 10 of Eichler book
# Quadratische Formen und orthogonale gruppen
# Also used the Scattone thesis as basis.
#
# Note: There are other works based on Eichler work that are analogous:
# ---Wall CTC, On The Orthogonal Groups of Unimodular Quadratic Forms
# ---James D.G. Indefinite Quadratic Forms of Determinant pm2p
# They also require the same hypothesis of having two hyperplanes.
INDEF_FORM_EichlerCriterion_TwoHyperplanesEven:=function(Qmat)
    local n, Block11, Block12, Gmat, GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative, EnumerateVectorOverDiscriminant, ListClasses, GRPstart, SetListClassesOrbitwise, GetApproximateGroup;
    n:=Length(Qmat);
    Block11:=SubBlock(Qmat, [1..4], [1..4]);
    Block12:=SubBlock(Qmat, [1..4], [5..n]);
    Gmat:=SubBlock(Qmat, [5..n], [5..n]);
    if Block11<>[[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]] or Block12<>NullMat(4,n-4) then
        Error("We do not have the two hyperplanes");
    fi;
    if INDEF_FORM_IsEven(Qmat)=false then
        Error("The lattice is not even");
    fi;
    # ComputeClasses
    ListClasses:=GetTranslationClasses(Gmat);
    GRPstart:=Group([IdentityMat(n)]);
    SetListClassesOrbitwise:=function(GRPmatr)
        local ListClassesExt, ListMatrGens, iMatr, GetPosition, ListPermGens, eMatrGen, eList, GRPperm, O, ListRepr;
        GRPstart:=GRPmatr;
        ListClassesExt:=List(ListClasses, x->Concatenation([0,0,0,0],x));
        GetPosition:=function(eVect)
            local i;
            for i in [1..Length(ListClassesExt)]
            do
                if SolutionIntMat(Qmat, ListClassesExt[i] - eVect)<>fail then
                    return i;
                fi;
            od;
            Error("Should not be there");
        end;
        ListPermGens:=[];
        ListMatrGens:=GeneratorsOfGroup(GRPmatr);
        for iMatr in [1..Length(ListMatrGens)]
        do
            Print("iMatr=", iMatr, " / ", Length(ListMatrGens), "\n");
            eMatrGen:=ListMatrGens[iMatr];
            eList:=List(ListClassesExt, x->GetPosition(x * eMatrGen));
            Add(ListPermGens, PermList(eList));
        od;
        GRPperm:=Group(ListPermGens);
        O:=Orbits(GRPperm, [1..Length(ListClasses)], OnPoints);
        Print("|O|=", List(O, Length), "\n");
        ListRepr:=List(O, x->x[1]);
        ListClasses:=ListClasses{ListRepr};
    end;
    GetApproximateGroup:=function()
        local GetLeftMultiplication, GetRightMultiplication, ExpandGenerator, ListGenerators, FuncInsert, eGenS, eGenT, eGen, eVect, i, eProd, BasisOrth, eOrth;
        # The quadratic form 2 x1 x2 + 2 x3 x4
        # We formulate it as (1/2) det U(x1,x2,x3,x4)
        # with U(x1,x2,x3,x4) = [ [x1 , -x4] , [x3 , x2] ]
        # If we multiply by a matrix in SL(2,Z) on the left and right
        GetLeftMultiplication:=function(P)
            local a, b, c, d, RetMat;
            # We compute the product [[a, b],[c,d]] U(x1,x2,x3,x4)
            # This gives us
            # [[a x1 + b x3 , -a x4 + b x2] , [ c x1 + d x3 , -c x4 + d x2]]
            # So Y = XA,
            # Y1 =  a x1 + b x3
            # Y2 = -c x4 + d x2
            # Y3 =  c x1 + d x3
            # Y4 =  a x4 - b x2
            a:=P[1][1];
            b:=P[1][2];
            c:=P[2][1];
            d:=P[2][2];
            RetMat:=[[a , 0, c, 0],
                     [0 , d, 0,-b],
                     [b , 0, d, 0],
                     [0 ,-c, 0, a]];
            return RetMat;
        end;
        GetRightMultiplication:=function(P)
            local a, b, c, d, RetMat;
            # We compute the product [ [x1 , -x4] , [x3 , x2] ]      [[a, b],[c,d]]
            # This gives us
            # [[a x1 - c x4 , b x1 - d x4] , [a x3 + c x2 , b x3 + d x2]]
            # So Y = XA,
            # Y1 = a x1 - c x4
            # Y2 = b x3 + d x2
            # Y3 = a x3 + c x2
            # Y4 = -b x1 + d x4
            a:=P[1][1];
            b:=P[1][2];
            c:=P[2][1];
            d:=P[2][2];
            RetMat:=[[ a, 0, 0,-b],
                     [ 0, d, c, 0],
                     [ 0, b, a, 0],
                     [-c, 0, 0, d]];
            return RetMat;
        end;
        ExpandGenerator:=function(eGen)
            local BigMat, i, j;
            BigMat:=IdentityMat(n);
            for i in [1..4]
            do
                for j in [1..4]
                do
                    BigMat[i][j]:=eGen[i][j];
                od;
            od;
            return BigMat;
        end;
        ListGenerators:=[];
        FuncInsert:=function(TheMat)
            if IsIntegralMat(TheMat)=false then
                Error("The matrix TheMat is not integral");
            fi;
            if TheMat * Qmat * TransposedMat(TheMat) <> Qmat then
                Error("The matrix is not preserving the quadratic form");
            fi;
            Add(ListGenerators, TheMat);
        end;
        for eGen in GeneratorsOfGroup(GRPstart)
        do
            FuncInsert(eGen);
        od;
        # Generators of SL2(Z)
        eGenS:=[[0,-1],[1,0]];
        eGenT:=[[1,1],[0,1]];
        for eGen in [eGenS, eGenT]
        do
            FuncInsert(ExpandGenerator(GetLeftMultiplication(eGen)));
            FuncInsert(ExpandGenerator(GetRightMultiplication(eGen)));
        od;
        # Now looking at the isotropic vectors, generating the Eichler transvections
        for i in [1..4]
        do
            eVect:=ListWithIdenticalEntries(n,0);
            eVect[i]:=1;
            eProd:=eVect * Qmat;
            BasisOrth:=NullspaceIntMat(TransposedMat([eProd]));
            for eOrth in BasisOrth
            do
                FuncInsert(INDEF_FORM_Eichler_Transvection(Qmat, eVect, eOrth));
            od;
        od;
        return Group(ListGenerators);
    end;
    # Notion of divisor
    # For an element y in a lattice L, define div(y) the integer such that
    # (L,y) = div(y) Z.
    # For each v in L* we can find smallest d such that w = dv in L.
    # And we then have div(w) = d
    EnumerateVectorOverDiscriminant:=function(X)
        local Ginv, Xdiv2, ListSolution, eClass1, eClass2, DivD, eClass3, X_res1, X_res2, Aclass3_norm, Aclass3_res1, Aclass3_res2, eClass4, Aclass4_norm, Aclass4_res2, eRec, quot, w, eSolution, diff, u, eNorm, DivD_frac, eProd;
        # For each vector v in M* / M
        # We want (x, v) divided by d.
        #
        # We should have e_I.x = d alpha_I
        # We should have (x,x) = X
        # For the first coordinates, we can write them as (0 , 0 , d , du)
        # By changing the value of u, we can shift the value of u which shifts the value
        # by 2 d^2, + or -.
        # The lattice is even so the 2 is not a problem.
        #
        # So, if we have some X, we want to find u\in Z and v in Z^{n-4} such that X = 2d^2 u + A[v]
        # with A the Gram matrix of M.
        # v must actually be of the form v = v0 + d A^{-1} h.
        # We are only interested on the value modulo 2d^2 of A on those vectors. So, it is a finite problem.
        # But it is a hard one as the search space might be very large.
        #
        if X mod 2 = 1 then
            return [];
        fi;
        Xdiv2 := X / 2;
        # The two hyperbolic planes are unimodular so for the discriminant, we only need
        # to consider the Gmat part.
        Ginv:=Inverse(Gmat);
        Print("Det(Gmat)=", DeterminantMat(Gmat), "\n");
        ListSolution:=[];
        for eClass1 in ListClasses
        do
            # The vector eClass1 is defined modulo an element of Gmat Z^n
            # Thus eClass2 is defined modulo an element of Z^n. This is an elements of L*
            # Thus eClass3 is defined modulo an element of d Z^n
            # So, we can write v = eClass3 + d w
            # which gets us A[v] = A[eClass3] + 2d w^T Gmat eClass3 + d^2 A[w]
            # So, the problem is actually linear.
            eClass2:=eClass1 * Ginv;
            DivD_frac:=RemoveFractionPlusCoef(eClass2).TheMult;
            DivD:=NumeratorRat(DivD_frac);
            eClass3:=eClass2 * DivD;
            #
            X_res1:=Xdiv2 mod DivD;
            X_res2:=Xdiv2 mod (DivD * DivD);
            Aclass3_norm:=(eClass3 * Gmat * eClass3) / 2;
            Aclass3_res1:=Aclass3_norm mod DivD;
            Aclass3_res2:=Aclass3_norm mod (DivD * DivD);
            if Aclass3_res1=X_res1 then # if not we cannot find a solution
                # and so
                # (X - A[eClass3])/2 = 2d w^T Gmat eClass3 + ....
                diff:=(X_res2 - Aclass3_res2) / DivD;
                eProd:=eClass3 * Gmat;
                eRec:=GcdVector(eProd);
                if eRec.TheGcd=0 then
                    quot:=0;
                else
                    quot:=diff / eRec.TheGcd;
                fi;
                if IsInt(quot) then
                    w:=quot * eRec.ListCoef;
                    eClass4:=eClass3 + DivD * w;
                    Aclass4_norm:=(eClass4 * Gmat * eClass4) / 2;
                    Aclass4_res2:=Aclass4_norm mod (DivD * DivD);
                    if Aclass4_res2<>X_res2 then
                        Error("A bug to resolve");
                    fi;
                    u:=(Xdiv2 - Aclass4_norm) / (DivD * DivD);
                    eSolution:=Concatenation([0,0,DivD, DivD*u], eClass4);
                    eNorm:=eSolution*Qmat*eSolution;
                    if eNorm<>X then
                        Error("A bug to resolve");
                    fi;
                    Add(ListSolution, eSolution);
                fi;
            fi;
        od;
        return ListSolution;
    end;
    GetCoveringOrbitRepresentatives:=function(X)
        local ListDiv, ListSolution, eDiv, eSol;
        if X = 0 then
            return EnumerateVectorOverDiscriminant(0);
        fi;
        ListDiv:=GetSquareDivisors(X);
        ListSolution:=[];
        for eDiv in ListDiv
        do
            for eSol in EnumerateVectorOverDiscriminant(X / (eDiv*eDiv))
            do
                Add(ListSolution, eDiv * eSol);
            od;
        od;
        return ListSolution;
    end;
    GetOneOrbitRepresentative:=function(X)
        local Xdiv;
        if X mod 2 = 1 then
            return fail;
        fi;
        Xdiv:=X / 2;
        return Concatenation([1,Xdiv], ListWithIdenticalEntries(n-2, 0));
    end;
    return rec(ApproximateGroup:=GetApproximateGroup(),
               SetListClassesOrbitwise:=SetListClassesOrbitwise,
               GetCoveringOrbitRepresentatives:=GetCoveringOrbitRepresentatives,
               GetOneOrbitRepresentative:=GetOneOrbitRepresentative);
end;


LORENTZ_ApproximateModel:=function(LorMat)
    local RecInput, ListFamily, GetApproximateGroup, SetListClassesOrbitwise, GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative;
    RecInput:=rec(TheOption:="total");
    ListFamily:=LORENTZ_EnumeratePerfect_DelaunayScheme(LorMat, RecInput);
    GetApproximateGroup:=function()
        return LORENTZ_OrthogonalGroup_Kernel(LorMat, ListFamily);
    end;
    SetListClassesOrbitwise:=function(GRPmatr)
        Error("Need to complete the code for SetListClassesOrbitwise");
    end;
    GetCoveringOrbitRepresentatives:=function(X)
        if X <> 0 then
            Error("X should be equal to 0. For others, please insert the code");
        fi;
        return LORENTZ_CandidateIsotropicVectors_Kernel(LorMat, ListFamily);
    end;
    GetOneOrbitRepresentative:=function(X)
        local ListCand;
        if X <> 0 then
            Error("X should be equal to 0. For others, please insert the code");
        fi;
        ListCand:=LORENTZ_CandidateIsotropicVectors_Kernel(LorMat, ListFamily);
        if Length(ListCand)=0 then
            return fail;
        fi;
        return ListCand[1];
    end;
    return rec(ApproximateGroup:=GetApproximateGroup(),
               SetListClassesOrbitwise:=SetListClassesOrbitwise,
               GetCoveringOrbitRepresentatives:=GetCoveringOrbitRepresentatives,
               GetOneOrbitRepresentative:=GetOneOrbitRepresentative);
end;


GetHyperbolicPlane:=function(Qmat)
    local n, ListVect, ThePerturb, HasSolution, MinScal, NbImprovement, ePerturb, M, eVect, NewVect, eScalPre, eScal, OneIso;
    n:=Length(Qmat);
    ListVect:=[INDEF_FindIsotropic(Qmat)];
#    Print("ListVect=", ListVect, "\n");
    ThePerturb:=IdentityMat(n);
    HasSolution:=false;
    MinScal:="unset";
    NbImprovement:=0;
    while(true)
    do
        ePerturb:=GetRandomMatrixPerturbation(n);
        ThePerturb:=ThePerturb * ePerturb;
        M:=ThePerturb * Qmat * TransposedMat(ThePerturb);
#        Print("Before INDEF_FindIsotropic\n");
#        PrintArray(M);
        eVect:=INDEF_FindIsotropic(M);
#        Print("After INDEF_FindIsotropic\n");
        NewVect:=eVect * ThePerturb;
        OneIso:=ListVect[1] * Inverse(ThePerturb);
#        Print("OneIso=", OneIso, "\n");
        if NewVect*Qmat*NewVect<>0 then
            Error("Vector should be isotrop");
        fi;
        if Position(ListVect, NewVect)=fail then
            for eVect in ListVect
            do
                eScalPre:=eVect * Qmat * NewVect;
                eScal:=AbsInt(eScalPre);
#                Print("eScal=", eScal, "\n");
                if eScal<>0 then
                    if HasSolution then
                        if eScal<MinScal then
                            MinScal:=eScal;
                            NbImprovement:=0;
                        else
                            if eScal=MinScal then
                                NbImprovement:=NbImprovement+1;
                                if NbImprovement=100 then
#                                    Print("Terminating MinScal=", MinScal, " scal=", eVect * Qmat * NewVect, "\n");
                                    if eScalPre>0 then
                                        return [eVect, NewVect];
                                    else
                                        return [eVect, -NewVect];
                                    fi;
                                fi;
                            fi;
                        fi;
                    else
                        HasSolution:=true;
                        MinScal:=eScal;
                    fi;
                fi;
            od;
            Add(ListVect, NewVect);
        fi;
    od;
end;


GetEichlerHyperplaneBasis:=function(Qmat)
    local Basis1, NSP, Qmat2, Basis2, HyperBasis, NSP2, FullBasis;
    Print("GetEichlerHyperplaneBasis, step 1\n");
    PrintArray(Qmat);
#    PrintArray(Qmat);
    Basis1:=GetHyperbolicPlane(Qmat);
#    Print("GetEichlerHyperplaneBasis, step 2\n");
    NSP:=NullspaceIntMat(TransposedMat(Basis1 * Qmat));
#    Print("GetEichlerHyperplaneBasis, step 3\n");
    Qmat2:=NSP * Qmat * TransposedMat(NSP);
#    Print("GetEichlerHyperplaneBasis, step 4\n");
#    PrintArray(Qmat2);
    Basis2:=GetHyperbolicPlane(Qmat2);
#    Print("GetEichlerHyperplaneBasis, step 5\n");
    HyperBasis:=Concatenation(Basis1, Basis2 * NSP);
#    Print("GetEichlerHyperplaneBasis, step 6\n");
    NSP2:=NullspaceIntMat(TransposedMat(HyperBasis * Qmat));
#    Print("GetEichlerHyperplaneBasis, step 7\n");
    FullBasis:=Concatenation(HyperBasis, NSP2);
    Print("GetEichlerHyperplaneBasis, step 8\n");
    return FullBasis;
end;


GetTwoPlaneEmbedding:=function(eBlock)
    local eEmbed, h1, h2, TwoPlanes, eProd;
    eEmbed:=IdentityMat(4);
    h1:=eBlock[1][2];
    h2:=eBlock[3][4];
    eEmbed[1][1]:=h1;
    eEmbed[3][3]:=h2;
    TwoPlanes:=[[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]];
    eProd:=eEmbed * TwoPlanes * TransposedMat(eEmbed);
    if eProd<>eBlock then
        Error("Failed to find the factorization");
    fi;
    return rec(eBlock:=TwoPlanes, eEmbed:=eEmbed);
end;


GetEasyIsometries:=function(Qmat)
    local n, GRA, i, j, LConn, ListQ, eQ, eConn, ListGenerators, f_insert, eConn1, eConn2, siz, GRPred, eGen, BigMat, u, v, dim, pos1, pos2, ListIdx, ListIdxNeg, len, LConnSel, idx;
    n:=Length(Qmat);
    GRA:=NullGraph(Group(()), n);
    for i in [1..n]
    do
        for j in [i+1..n]
        do
            if Qmat[i][j] <> 0 then
                AddEdgeOrbit(GRA, [i,j]);
                AddEdgeOrbit(GRA, [j,i]);
            fi;
        od;
    od;
    LConn:=ConnectedComponents(GRA);
    ListQ:=[];
    for eConn in LConn
    do
        eQ:=SubBlock(Qmat, eConn, eConn);
        Add(ListQ, eQ);
    od;
    ListGenerators:=[IdentityMat(n)];
    f_insert:=function(eP)
        if IsIntegralMat(eP)=false then
            Error("eP should be integral");
        fi;
        if eP * Qmat * TransposedMat(eP)<>Qmat then
            Error("eP does not preserve the quadratic form");
        fi;
        Add(ListGenerators, eP);
    end;
    # Insert the isometries of each block
    for i in [1..Length(LConn)]
    do
        eQ:=ListQ[i];
        eConn:=LConn[i];
        siz:=Length(eConn);
        if IsPositiveDefiniteSymmetricMatrix(eQ) then
            GRPred:=ArithmeticAutomorphismGroup([eQ]);
            for eGen in GeneratorsOfGroup(GRPred)
            do
                BigMat:=IdentityMat(n);
                for u in [1..siz]
                do
                    for v in [1..siz]
                    do
                        BigMat[eConn[u]][eConn[v]] := eGen[u][v];
                    od;
                od;
                f_insert(BigMat);
            od;
        fi;
    od;
    # Insert the permutations of each block
    for eQ in Set(ListQ)
    do
        dim:=Length(eQ);
        ListIdx:=Filtered([1..Length(ListQ)], x->ListQ[x] = eQ);
        ListIdxNeg:=Filtered([1..Length(ListQ)], x->ListQ[x] <> eQ);
        len:=Length(ListIdx);
        LConnSel:=LConn{ListIdx};
        for eGen in GeneratorsOfGroup(SymmetricGroup([1..len]))
        do
            BigMat:=NullMat(n,n);
            for i in [1..len]
            do
                j:=OnPoints(i,eGen);
                eConn1:=LConnSel[i];
                eConn2:=LConnSel[j];
                for u in [1..dim]
                do
                    pos1:=eConn1[u];
                    pos2:=eConn2[u];
                    BigMat[pos1][pos2]:=1;
                od;
            od;
            for idx in ListIdxNeg
            do
                for u in LConn[idx]
                do
                    BigMat[u][u]:=1;
                od;
            od;
            f_insert(BigMat);
        od;
    od;
    return Group(ListGenerators);
end;


INDEF_FORM_GetApproximateModel:=function(Qmat)
    local n, FullBasis, QmatRed, Block11, Block12, Block22, TwoPlanes, QmatExt, RecApprox, eEmbed, ListCoset, RecStab_RightCoset, TheRec, ListGenerators, eGen, fGen, NewGen, GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative, Hmat, QmatEichler, TheEmbed, GRPeasy, TheRecEmbed, eProdEmbed;
    if LORENTZ_IsLorentzian(Qmat) then
        return LORENTZ_ApproximateModel(Qmat);
    fi;
    # We try the Eichler model as it is the thing that works
    n:=Length(Qmat);
    Print("Before GetEichlerHyperplaneBasis\n");
    FullBasis:=GetEichlerHyperplaneBasis(Qmat);
    Print("After GetEichlerHyperplaneBasis\n");
    if AbsInt(DeterminantMat(FullBasis))<>1 then
        Print("Further works is needed here");
    fi;
    QmatRed:=FullBasis * Qmat * TransposedMat(FullBasis);
    Block11:=SubBlock(QmatRed, [1..4], [1..4]);
    Block12:=SubBlock(QmatRed, [1..4], [5..n]);
    Block22:=SubBlock(QmatRed, [5..n], [5..n]);
    TwoPlanes:=[[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]];
    if Block11=TwoPlanes then
        Print("INDEF_FORM_GetApproximateModel, case 1\n");
        TheRec:=INDEF_FORM_EichlerCriterion_TwoHyperplanesEven(QmatRed);
        GRPeasy:=GetEasyIsometries(QmatRed);
        TheRec.SetListClassesOrbitwise(GRPeasy);
        ListGenerators:=[];
        for eGen in GeneratorsOfGroup(TheRec.ApproximateGroup)
        do
            NewGen:=Inverse(FullBasis) * eGen * FullBasis;
            Add(ListGenerators, NewGen);
        od;
        GetCoveringOrbitRepresentatives:=function(X)
            return List(TheRec.GetCoveringOrbitRepresentatives(X), x->x*FullBasis);
        end;
        GetOneOrbitRepresentative:=function(X)
            local eAns;
            eAns:=TheRec.GetOneOrbitRepresentative(X);
            if eAns=fail then
                return fail;
            fi;
            return eAns * FullBasis;
        end;
        return rec(ApproximateGroup:=Group(ListGenerators),
                   GetCoveringOrbitRepresentatives:=GetCoveringOrbitRepresentatives,
                   GetOneOrbitRepresentative:=GetOneOrbitRepresentative);
    fi;
    TheRecEmbed:=GetTwoPlaneEmbedding(Block11);
    if TheRecEmbed.eBlock=TwoPlanes and Block12=NullMat(4,n-4) then
        Print("INDEF_FORM_GetApproximateModel, case 2\n");
        QmatExt:=LORENTZ_AssembleDiagBlocks([TwoPlanes, Block22]);
        Print("QmatExt=", QmatExt, "\n");
        eEmbed:=LORENTZ_AssembleDiagBlocks([TheRecEmbed.eEmbed, IdentityMat(n-4)]);
        if eEmbed * QmatExt * TransposedMat(eEmbed) <> QmatRed then
            Error("eEmned fails");
        fi;
        eProdEmbed:=Inverse(FullBasis) * eEmbed;
        if eProdEmbed * QmatExt * TransposedMat(eProdEmbed) <> Qmat then
            Error("eProdEmbed fails");
        fi;
        Print("eEmbed=", eEmbed, "\n");
        Print("Before INDEF_FORM_GetApproximateModel, case 1\n");
        RecApprox:=INDEF_FORM_GetApproximateModel(QmatExt);
        Print("We have RecApprox\n");
        RecStab_RightCoset:=LinearSpace_Stabilizer_RightCoset(RecApprox.ApproximateGroup, eEmbed);
        Print("We have RecStab_RightCoset\n");
        ListCoset:=LinearSpace_ExpandListListCoset(n, RecStab_RightCoset.ListListCoset);
        Print("We have ListCoset\n");
        ListGenerators:=[IdentityMat(n)];
        for eGen in GeneratorsOfGroup(RecStab_RightCoset.GRPmatr)
        do
            fGen:=eProdEmbed * eGen * Inverse(eProdEmbed);
            if IsIntegralMat(fGen)=false then
                Error("fGen should be integral");
            fi;
            if fGen * Qmat * TransposedMat(fGen)<>Qmat then
                Error("fGen does not preserve Qmat scalar product");
            fi;
            Add(ListGenerators, fGen);
        od;
        Print("We have ListGenerators\n");
        GetCoveringOrbitRepresentatives:=function(X)
            local ListRepr, eRepr, eCos, fRepr, eSol, fSol;
            ListRepr:=[];
            for eRepr in RecApprox.GetCoveringOrbitRepresentatives(X)
            do
                for eCos in ListCoset
                do
                    fRepr:=eRepr * eCos;
                    eSol:=SolutionIntMat(eEmbed, fRepr);
                    if eSol<>fail then
                        fSol:=eSol * FullBasis;
                        if fSol * Qmat * fSol <> X then
                            Error("The vector is not of norm X");
                        fi;
                        Add(ListRepr, fSol);
                    fi;
                od;
            od;
            return ListRepr;
        end;
        GetOneOrbitRepresentative:=function(X)
            local eRepr, eCos, fRepr, eSol, fSol, ListRepr;
            eRepr:=RecApprox.GetOneOrbitRepresentative(X);
            if eRepr=fail then
                return fail;
            fi;
            for eCos in ListCoset
            do
                fRepr:=eRepr * eCos;
                eSol:=SolutionIntMat(eEmbed, fRepr);
                if eSol<>fail then
                    fSol:=eSol * FullBasis;
                    if fSol * Qmat * fSol <> X then
                        Error("The vector is not of norm X");
                    fi;
                    return fSol;
                fi;
            od;
            # We hoped to avoid that call, but here goes
            ListRepr:=GetCoveringOrbitRepresentatives(X);
            if Length(ListRepr) > 0 then
                return ListRepr[1];
            fi;
            return fail;
        end;
        return rec(ApproximateGroup:=Group(ListGenerators),
                   GetCoveringOrbitRepresentatives:=GetCoveringOrbitRepresentatives,
                   GetOneOrbitRepresentative:=GetOneOrbitRepresentative);
    fi;
    Error("Something more to write");
end;


ExpandMatrix:=function(eEndoRed)
    local n, TheBigMat, i, j;
    n:=Length(eEndoRed) + 1;
    TheBigMat:=IdentityMat(n);
    for i in [1..n-1]
    do
        for j in [1..n-1]
        do
            TheBigMat[i][j]:=eEndoRed[i][j];
        od;
    od;
    return TheBigMat;
end;


# Isometry group defined on a p dimensional space for a quadratic form Qp.
# We extend the quadratic form to dimension n with
# Qn = | Qp 0 |
#      | 0  0 |
# If the original matrices satisfy Pp Qp Pp^T = Qp
# the the extended matrices must satisfy Pn Qn Pn^T = Qn
# So, in block formulation
# Pn = | A B |
#      | C D |
# and so
# Pn Qn Pn^T = | A B |     | Qp 0 |     | A^T C^T |
#              | C D |  x  | 0  0 |  x  | B^T D^T |
#            = | A Qp 0 |     | A^T C^T |
#              | C Qp 0 |  x  | B^T D^T |
#            = | A Qp A^T  A Qp C^T |
#              | C Qp A^T  C Qp C^T |
#            = | Qp 0 |
#              | 0  0 |
# And so we get A Qp A^T = Qp , C Qp A^T = 0 , C Qp C^T = 0
# The equation A Qp A^T forces A to be an isometry of Qp.
# The equation C Qp A^T forces C to be 0 from which the rest follows.
ExtendIsometryGroup:=function(GRPmatr, p, n)
    local ListGens, eGen, NewMat, i, j;
    ListGens:=[];
    for eGen in GeneratorsOfGroup(GRPmatr)
    do
        NewMat:=IdentityMat(n);
        for i in [1..p]
        do
            for j in [1..p]
            do
                NewMat[i][j]:=eGen[i][j];
            od;
        od;
        Add(ListGens, NewMat);
    od;
    if n > p then
        for eGen in GeneratorsOfGroup(GeneralLinearGroup(n-p,Integers))
        do
            NewMat:=IdentityMat(n);
            for i in [1..n-p]
            do
                for j in [1..n-p]
                do
                    NewMat[i+p][j+p]:=eGen[i][j];
                od;
            od;
            Add(ListGens, NewMat);
        od;
        for i in [1..p]
        do
            NewMat:=IdentityMat(n);
            NewMat[i][p+1]:=1;
            Add(ListGens, NewMat);
        od;
    fi;
    return Group(ListGens);
end;


# Suppose we have a flag formed by vectors f1 = e1, f2 = {e1, e2}, f3 = {e1, e2, e3}
# We want to find a generating set of the automorphism group
GetAutomorphismOfFlag:=function(n)
    local LGen, TheMat, i, j;
    LGen:=[];
    for i in [1..n]
    do
        TheMat:=IdentityMat(n);
        TheMat[i][i]:=-1;
        Add(LGen, TheMat);
    od;
    for i in [1..n]
    do
        for j in [1..i-1]
        do
            TheMat:=IdentityMat(n);
            TheMat[i][j]:=1;
            Add(LGen, TheMat);
        od;
    od;
    return Group(LGen);
end;


# It is the same as above but the group in question
# is only the flag presrving one
ExtendIsometryGroup_Triangular:=function(GRPmatr, p, n)
    local ListGens, eGen, NewMat, i, j, SubNPgroup, idx;
    ListGens:=[];
    for eGen in GeneratorsOfGroup(GRPmatr)
    do
        NewMat:=IdentityMat(n);
        for i in [1..p]
        do
            for j in [1..p]
            do
                NewMat[i][j]:=eGen[i][j];
            od;
        od;
        Add(ListGens, NewMat);
    od;
    SubNPgroup:=GetAutomorphismOfFlag(n-p);
    for eGen in GeneratorsOfGroup(SubNPgroup)
    do
        NewMat:=IdentityMat(n);
        for i in [1..n-p]
        do
            for j in [1..n-p]
            do
                NewMat[i+p][j+p]:=eGen[i][j];
            od;
        od;
        Add(ListGens, NewMat);
    od;
    for i in [1..p]
    do
        for idx in [p+1..n]
        do
            NewMat:=IdentityMat(n);
            NewMat[i][idx]:=1;
            Add(ListGens, NewMat);
        od;
    od;
    return Group(ListGens);
end;


f_get_list_spaces:=function(ListVect)
    local len, ListSpaces, i;
    len:=Length(ListVect);
    ListSpaces:=[];
    for i in [1..len]
    do
        Add(ListSpaces, ListVect{[1..i]});
    od;
    return ListSpaces;
end;


INDEF_FORM_GetVectorStructure:=function(Qmat, v)
    local n, eNorm, eProd, NSP, vDisc, Pmat, PmatInv, i, j, GramMatRed, MapOrthogonalSublatticeGroup, MapOrthogonalSublatticeEndomorphism;
    n:=Length(Qmat);
    eNorm:=v * Qmat * v;
    eProd:=v*Qmat;
    NSP:=NullspaceIntMat(TransposedMat([eProd]));
    GramMatRed:=NSP*Qmat*TransposedMat(NSP);
    Pmat:="unset";
    PmatInv:="unset";
    if eNorm<>0 then
        Pmat:=NullMat(n,n);
        for i in [1..n-1]
        do
            for j in [1..n]
            do
                Pmat[i][j]:=NSP[i][j];
            od;
        od;
        for j in [1..n]
        do
            Pmat[n][j]:=v[j];
        od;
        Print("INDEF_FORM_GetVectorStructure |Pmat|=", DeterminantMat(Pmat), "\n");
        PmatInv:=Inverse(Pmat);
    fi;
    MapOrthogonalSublatticeEndomorphism:=function(eEndoRed)
        local TheBigMat, RetMat, Subspace1, Subspace2, vImg;
        if eNorm<>0 then
            TheBigMat:=ExpandMatrix(eEndoRed);
            RetMat:=PmatInv * TheBigMat * Pmat;
            if v * RetMat <> v then
                Error("Vector should be invariant (general case)");
            fi;
        else
            Subspace1:=eEndoRed * NSP;
            Subspace2:=NSP;
            Print("Before calling LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1 in INDEF_FORM_GetVectorStructure\n");
            RetMat:=LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(Qmat, Subspace1, Qmat, Subspace2);
            vImg:=v * RetMat;
            if vImg <> v and vImg <> -v then
                Error("Vector should be invariant (isotropic case)");
            fi;
            if vImg = -v then
                RetMat:=-RetMat;
            fi;
        fi;
        return RetMat;
    end;
    MapOrthogonalSublatticeGroup:=function(GRPmatr)
        local GRPret;
        GRPret:=PersoGroupMatrix(List(GeneratorsOfGroup(GRPmatr), MapOrthogonalSublatticeEndomorphism), n);
        return GRPret;
    end;
    return rec(MapOrthogonalSublatticeEndomorphism:=MapOrthogonalSublatticeEndomorphism,
               MapOrthogonalSublatticeGroup:=MapOrthogonalSublatticeGroup,
               Pmat:=Pmat,
               PmatInv:=PmatInv,
               GramMatRed:=GramMatRed,
               NSP:=NSP);
end;


#
#
# The automorphism on the Plane^{perp} can be extended to full autmorphism.
# But if dim(Plane) > 1 then the extension is not unique.
# What we need is to find some lattice such that any extension
#
INDEF_FORM_GetRec_IsotropicKplane:=function(Qmat, Plane)
    local dimSpace, dim, eProd, WitnessIsotropy, NSP, eLine, eSol, GramMatRed, MapOrthogonalSublatticeEndomorphism, MapOrthogonalSublatticeGroup, check_generator;
    dimSpace:=Length(Qmat);
    dim:=Length(Plane);
    eProd:=Plane * Qmat;
    WitnessIsotropy:=Plane * Qmat * TransposedMat(Plane);
    if WitnessIsotropy<>NullMat(dim,dim) then
        Error("The matrix Plane does not define a totally isotropic space");
    fi;
    NSP:=NullspaceIntMat(TransposedMat(eProd));
    GramMatRed:=NSP * Qmat * TransposedMat(NSP);
    check_generator:=function(eEndoRed, RetMat)
        local PlaneImg, TransRed, RetMat_red;
        PlaneImg:=Plane * RetMat;
        TransRed:=List(PlaneImg, x->SolutionMat(Plane, x));
        if AbsInt(DeterminantMat(TransRed))<>1 then
            Error("TransRed should have absolute determinant 1\n");
        fi;
#        Print("  TransRed=", TransRed, "\n");
        if TestEqualitySpace(PlaneImg, Plane)=false then
            Error("Plane should be invariant (isotropic case)");
        fi;
        RetMat_red:=List(NSP, x->SolutionMat(NSP, x * RetMat));
        if RetMat_red<>eEndoRed then
            Error("RetMat_red restricted to the nullspace is not the original eEndoRed");
        fi;
    end;
    #
    # We map automorphism group of a sublattice to the full group.
    # The difficulty os that the rational kernel is of strictly positive Q-rank
    # if dim(Plane) > 1.
    MapOrthogonalSublatticeGroup:=function(GRPmatr)
        local Subspace1, ListGenTotal, ListD, eEndoRed, Subspace2, TheRec, RetMat, TheDen, TheKer;
        Subspace1:=NSP;
        ListGenTotal:=[];
        ListD:=[ 1 ];
        for eEndoRed in GeneratorsOfGroup(GRPmatr)
        do
            Subspace2:=eEndoRed * NSP;
            TheRec:=LORENTZ_ExtendOrthogonalIsotropicIsomorphism(Qmat, Subspace1, Qmat, Subspace2);
            RetMat:=TheRec.get_one_transformation();
            check_generator(eEndoRed, RetMat);
            Add(ListGenTotal, RetMat);
            Add(ListD, GetDenominatorMatrix(RetMat));
        od;
        TheDen:=Lcm(ListD);
        Print("ListD=", ListD, " TheDen=", TheDen, "\n");
        TheRec:=LORENTZ_ExtendOrthogonalIsotropicIsomorphism(Qmat, Subspace1, Qmat, Subspace1);
        TheKer:=TheRec.get_kernel_generating_set(TheDen);
#        Print("TheKer=", TheKer, "\n");
        eEndoRed:=IdentityMat(Length(NSP));
        for RetMat in GeneratorsOfGroup(TheKer)
        do
            check_generator(eEndoRed, RetMat);
            Add(ListGenTotal, RetMat);
        od;
        return PersoGroupMatrix(ListGenTotal, dimSpace);
    end;
    return rec(Plane:=Plane, Qmat:=Qmat,
               MapOrthogonalSublatticeGroup:=MapOrthogonalSublatticeGroup,
               GramMatRed:=GramMatRed,
               NSP:=NSP);
end;


GetFirstNorm:=function(ApproxModel)
    local X, eVect;
    X:=1;
    while(true)
    do
        eVect:=ApproxModel.GetOneOrbitRepresentative(X);
        if eVect<>fail then
            return rec(X:=X, eVect:=eVect);
        fi;
        X:=X+1;
    od;
end;


IndefiniteReductionTrivial:=function(M)
    return rec(Mred:=M, B:=IdentityMat(Length(M)));
end;


# Let us think about computing with some stabilizer.
#
INDEF_FORM_Machinery_AllFct:=function()
    local ListResultsIsomorphism, GetListVert, GetGraListIConn, GetEquivalence, DatabaseKeyValue,
          INDEF_FORM_GetOrbitRepresentative, INDEF_FORM_TestEquivalence, INDEF_FORM_EquivalenceVector,
          INDEF_FORM_AutomorphismGroup, INDEF_FORM_StabilizerVector,
          INDEF_FORM_AutomorphismGroup_Memoize, INDEF_FORM_TestEquivalence_Memoize,
          INDEF_FORM_AutomorphismGroup_Reduction, INDEF_FORM_TestEquivalence_Reduction,
          INDEF_FORM_AutomorphismGroup_Kernel, INDEF_FORM_TestEquivalence_Kernel,
          Kernel_Equivalence_Qmat,
          INDEF_FORM_Equivalence_IsotropicKplane, INDEF_FORM_RightCosets_IsotropicKplane,
          INDEF_FORM_Stabilizer_IsotropicKplane, INDEF_FORM_GetOrbit_IsotropicKplane,
          INDEF_FORM_Equivalence_IsotropicKstuff_Kernel, INDEF_FORM_RightCosets_IsotropicKstuff_Kernel,
          INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel, INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel,
          INDEF_FORM_Invariant_IsotropicKstuff_Kernel, INDEF_FORM_Invariant_IsotropicKplane,
          f_equiv_plane, f_stab_plane, f_equiv_flag, f_stab_flag,
          INDEF_FORM_Equivalence_IsotropicKflag, INDEF_FORM_Invariant_IsotropicKflag,
          INDEF_FORM_Stabilizer_IsotropicKflag, INDEF_FORM_RightCosets_IsotropicKflag,
          INDEF_FORM_GetOrbit_IsotropicKflag;
    # We store some classes of isomorphism, that is what we have computed.
    # Every time we compute some non-trivial isomorphism, we store it in the list of results.
    # The list is not refined for isomorphisms because curating the database is itself
    # an additional work.
    ListResultsIsomorphism:=[];
    # Assumes M1 and M2 are in the graph
    GetListVert:=function()
        local ListVertM, eRec;
        ListVertM:=[];
        for eRec in ListResultsIsomorphism
        do
            AddSet(ListVertM, eRec.M1);
            AddSet(ListVertM, eRec.M2);
        od;
        return ListVertM;
    end;
    GetGraListIConn:=function(ListVertM)
        local GRA, eRec, posA, posB, LConn, ListIConn, iConn, eConn;
        GRA:=NullGraph(Group(()), Length(ListVertM));
        for eRec in ListResultsIsomorphism
        do
            if eRec.eEquiv<>fail then
                posA:=Position(ListVertM, eRec.M1);
                posB:=Position(ListVertM, eRec.M2);
                AddEdgeOrbit(GRA, [posA, posB]);
                AddEdgeOrbit(GRA, [posB, posA]);
            fi;
        od;
        LConn:=ConnectedComponents(GRA);
        ListIConn:=ListWithIdenticalEntries(Length(ListVertM),0);
        for iConn in [1..Length(LConn)]
        do
            eConn:=LConn[iConn];
            ListIConn{eConn}:=ListWithIdenticalEntries(Length(eConn), iConn);
        od;
        return rec(GRA:=GRA, ListIConn:=ListIConn);
    end;
    GetEquivalence:=function(GRA, ListVertM, M1, M2)
        local eP, pos1, pos2, TheShort, i, posA, posB, FindEquiv, eRec, posX, posY;
        eP:=IdentityMat(Length(M1));
        pos1:=Position(ListVertM, M1);
        pos2:=Position(ListVertM, M2);
        TheShort:=FindShortestPath(GRA, pos1, pos2);
        for i in [2..Length(TheShort)]
        do
            posA:=TheShort[i-1];
            posB:=TheShort[i];
            FindEquiv:=false;
            for eRec in ListResultsIsomorphism
            do
                posX:=Position(ListVertM, eRec.M1);
                posY:=Position(ListVertM, eRec.M2);
                if posX=posA and posY=posB then
                    eP:=eRec.eEquiv * eP;
                    FindEquiv:=true;
                fi;
                if posX=posB and posY=posA then
                    eP:=Inverse(eRec.eEquiv) * eP;
                    FindEquiv:=true;
                fi;
            od;
            if FindEquiv=false then
                Error("We failed to find a locl equivalence");
            fi;
        od;
        if eP * M1 * TransposedMat(eP) <> M2 then
            Error("The construction of equivalence has failed");
        fi;
        return eP;
    end;
    INDEF_FORM_TestEquivalence_Memoize:=function(M1, M2)
        local ListVertM, pos1, pos2, ePairGra, iConn1, iConn2, eRec, posA, posB, iConnA, iConnB, eEquiv;
        if M1 = M2 then
            return IdentityMat(Length(M1));
        fi;
        ListVertM:=GetListVert();
        pos1:=Position(ListVertM, M1);
        pos2:=Position(ListVertM, M2);
        if pos1<>fail and pos2<>fail then
            ePairGra:=GetGraListIConn(ListVertM);
            iConn1:=ePairGra.ListIConn[pos1];
            iConn2:=ePairGra.ListIConn[pos2];
            if iConn1<>iConn2 then
                for eRec in ListResultsIsomorphism
                do
                    if eRec.eEquiv=fail then
                        posA:=Position(ListVertM, eRec.M1);
                        posB:=Position(ListVertM, eRec.M2);
                        iConnA:=ePairGra.ListIConn[posA];
                        iConnB:=ePairGra.ListIConn[posB];
                        # We have already one non-isomorphism. Returns fail
                        if Set([iConn1,iConn2]) = Set([iConnA,iConnB]) then
                            return fail;
                        fi;
                    fi;
                od;
                # We could not conclude tha M1 and M2 are not isomorphic.
                # The next computation will decide of that
            else
                return GetEquivalence(ePairGra.GRA, ListVertM, M1, M2);
            fi;
        fi;
        eEquiv:=INDEF_FORM_TestEquivalence_Kernel(M1,M2);
        Add(ListResultsIsomorphism, rec(M1:=M1, M2:=M2, eEquiv:=eEquiv));
        return eEquiv;
    end;
    DatabaseKeyValue:=rec(ListKey:=[], ListValue:=[]);
    INDEF_FORM_AutomorphismGroup_Memoize:=function(M)
        local ListVertM, pos2, ePairGra, n_keyval, iKey, M1, pos1, GRP1, ListGens, eEquiv, GRP, eGen1, eGen;
        ListVertM:=GetListVert();
        pos2:=Position(ListVertM, M);
        if pos2<>fail then
            ePairGra:=GetGraListIConn(ListVertM);
            n_keyval:=Length(DatabaseKeyValue.ListKey);
            for iKey in [1..n_keyval]
            do
                M1:=DatabaseKeyValue.ListKey[iKey];
                pos1:=Position(ListVertM, M1);
                if pos1<>fail then
                    if ePairGra.ListIConn[pos1] = ePairGra.ListIConn[pos2] then
                        GRP1:=DatabaseKeyValue.ListValue[iKey];
                        eEquiv:=GetEquivalence(ePairGra.GRA, ListVertM, M1, M);
                        ListGens:=[];
                        for eGen1 in GeneratorsOfGroup(GRP1)
                        do
                            if eGen1 * M1 * TransposedMat(eGen1) <> M1 then
                                Error("eGen1 should preserve M1");
                            fi;
                            eGen:=eEquiv * eGen1 * Inverse(eEquiv);
                            if eGen * M * TransposedMat(eGen) <> M then
                                Error("eGen should preserve M");
                            fi;
                        od;
                        return PersoGroupMatrix(ListGens, Length(M));
                    fi;
                fi;
            od;
        fi;
        GRP:=INDEF_FORM_AutomorphismGroup_Kernel(M);
        Add(DatabaseKeyValue.ListKey, M);
        Add(DatabaseKeyValue.ListValue, GRP);
        return GRP;
    end;
    INDEF_FORM_TestEquivalence_Kernel:=function(Qmat1, Qmat2)
        local Block1, Block2, ApproxModel1, eRec1, X, v1, ApproxModel2, v2, test, ListCand;
        Block1:=INDEF_FORM_GetAttackScheme(Qmat1);
        Block2:=INDEF_FORM_GetAttackScheme(Qmat2);
        Print("Beginning of INDEF_FORM_TestEquivalence Block1.h=", Block1.h, " Block2.h=", Block2.h, "\n");
        if Block1.h = 0 then
            return INDEF_FORM_TestEquivalence_PosNeg(Qmat1, Qmat2);
        fi;
        if Block1.h = 1 then
            test:=LORENTZ_TestIsomorphismLorentzianMatrices(Block1.mat, Block2.mat);
            if test=false then
                return fail;
            fi;
            return test;
        fi;
        Print("Before INDEF_FORM_GetApproximateModel, case 2\n");
        ApproxModel1:=INDEF_FORM_GetApproximateModel(Qmat1);
        eRec1:=GetFirstNorm(ApproxModel1);
        X:=eRec1.X;
        v1:=eRec1.eVect;
        Print("X=", X, " v1=", v1, "\n");
        Print("Before INDEF_FORM_GetApproximateModel, case 3\n");
        ApproxModel2:=INDEF_FORM_GetApproximateModel(Qmat2);
        ListCand:=ApproxModel2.GetCoveringOrbitRepresentatives(X);
        for v2 in ListCand
        do
            test:=INDEF_FORM_EquivalenceVector(Qmat1, Qmat2, v1, v2);
            if test<>fail then
                return test;
            fi;
        od;
        Print("Found not isomorphic because of no vector equivalence found |ListCand|=", Length(ListCand), "\n");
        return fail;
    end;
    INDEF_FORM_TestEquivalence_Reduction:=function(Qmat1, Qmat2)
        local RecRed1, RecRed2, test, testNew;
        if RankMat(Qmat1)<>Length(Qmat1) or RankMat(Qmat2)<>Length(Qmat2) then
            Error("We should have Qmat1 and Qmat2 of full rank");
        fi;
        if INDEF_FORM_Invariant(Qmat1)<>INDEF_FORM_Invariant(Qmat2) then
            Print("Found not isomorphic because of different invariants\n");
            return fail;
        fi;
        # RecRed1:=IndefiniteReductionTrivial(Qmat1);
        # RecRed2:=IndefiniteReductionTrivial(Qmat2);
        RecRed1:=IndefiniteReduction(Qmat1);
        RecRed2:=IndefiniteReduction(Qmat2);
        # RecRed1.B * Qmat1 * TransposedMat(RecRed1.B) = RecRed1.Mred
        # RecRed2.B * Qmat2 * TransposedMat(RecRed2.B) = RecRed2.Mred
        test:=INDEF_FORM_TestEquivalence_Memoize(RecRed1.Mred, RecRed2.Mred);
        if test=fail then
            return fail;
        fi;
        # test * RedRed1.Mred * TransposedMat(test) = RedRed2.Mred
        # test * RecRed1.B * Qmat1 * TransposedMat(RecRed1.B) * TransposedMat(test) = RecRed2.B * Qmat2 * TransposedMat(RecRed2.B)
        testNew:=Inverse(RecRed2.B) * test * RecRed1.B;
        if testNew * Qmat1 * TransposedMat(testNew) <> Qmat2 then
            Error("We failed to have an equivalence");
        fi;
        return testNew;
    end;
    INDEF_FORM_TestEquivalence:=function(Qmat1, Qmat2)
        local NSP1, TheCompl1, FullBasis1, NSP2, TheCompl2, FullBasis2, n, p, i, j, QmatRed1, QmatRed2, test, TheEquivTest, TheEquiv;
        n:=Length(Qmat1);
        NSP1:=NullspaceIntMat(Qmat1);
        TheCompl1:=SubspaceCompletionInt(NSP1, n);
        FullBasis1:=Concatenation(TheCompl1, NSP1);
        NSP2:=NullspaceIntMat(Qmat2);
        TheCompl2:=SubspaceCompletionInt(NSP2, n);
        FullBasis2:=Concatenation(TheCompl2, NSP2);
        p:=Length(TheCompl1);
        QmatRed1:=TheCompl1 * Qmat1 * TransposedMat(TheCompl1);
        QmatRed2:=TheCompl2 * Qmat2 * TransposedMat(TheCompl2);
        test:=INDEF_FORM_TestEquivalence_Reduction(QmatRed1, QmatRed2);
        if test=fail then
            return fail;
        fi;
        TheEquivTest:=IdentityMat(n);
        for i in [1..p]
        do
            for j in [1..p]
            do
                TheEquivTest[i][j]:=test[i][j];
            od;
        od;
        TheEquiv:=Inverse(FullBasis2) * TheEquivTest * FullBasis1;
        if TheEquiv * Qmat1 * TransposedMat(TheEquiv)<>Qmat2 then
            Error("This is not an equivalence");
        fi;
        return TheEquiv;
    end;
    #
    # The automorphism group business
    #
    INDEF_FORM_AutomorphismGroup_Kernel:=function(Qmat)
        local eBlock, ListGenerators, ApproxModel, eRec, X, v1, v2, test, eGen, f_insert;
        eBlock:=INDEF_FORM_GetAttackScheme(Qmat);
        Print("Beginning of INDEF_FORM_AutomorphismGroup Block.h=", eBlock.h, "\n");
        if eBlock.h = 0 then
            return INDEF_FORM_AutomorphismGroup_PosNeg(Qmat);
        fi;
        if eBlock.h = 1 then
            return LORENTZ_OrthogonalGroup(eBlock.mat);
        fi;
        ListGenerators:=[];
        f_insert:=function(eGen)
            if eGen*Qmat*TransposedMat(eGen)<>Qmat then
                Error("The matrix generator does not preserve Qmat");
            fi;
            AddSet(ListGenerators, eGen);
        end;
        Print("Before INDEF_FORM_GetApproximateModel, case 4\n");
        ApproxModel:=INDEF_FORM_GetApproximateModel(Qmat);
        eRec:=GetFirstNorm(ApproxModel);
        X:=eRec.X;
        v1:=eRec.eVect;
        Print("X=", X, " v1=", v1, "\n");
        for eGen in GeneratorsOfGroup(ApproxModel.ApproximateGroup)
        do
            f_insert(eGen);
        od;
        for eGen in GeneratorsOfGroup(INDEF_FORM_StabilizerVector(Qmat, eRec.eVect))
        do
            f_insert(eGen);
        od;
        for v2 in ApproxModel.GetCoveringOrbitRepresentatives(X)
        do
            test:=INDEF_FORM_EquivalenceVector(Qmat, Qmat, v1, v2);
            if test<>fail then
                f_insert(eGen);
            fi;
        od;
        return PersoGroupMatrix(ListGenerators, Length(Qmat));
    end;
    INDEF_FORM_AutomorphismGroup_Reduction:=function(Qmat)
        local n, RecRed, GRP1, LGen2, eGen1, eGen2;
        if RankMat(Qmat)<>Length(Qmat) then
            Error("We should have Qmat of full rank");
        fi;
        # RecRed:=IndefiniteReductionTrivial(Qmat);
        RecRed:=IndefiniteReduction(Qmat);
        n:=Length(Qmat);
        # RecRed.B * Qmat1 * TransposedMat(RecRed.B) = RecRed.Mred
        GRP1:=INDEF_FORM_AutomorphismGroup_Memoize(RecRed.Mred);
        LGen2:=[];
        for eGen1 in GeneratorsOfGroup(GRP1)
        do
            eGen2:=Inverse(RecRed.B) * eGen1 * RecRed.B;
            if eGen2 * Qmat * TransposedMat(eGen2)<>Qmat then
                Error("The matrix is not an equivalence");
            fi;
            Add(LGen2, eGen2);
        od;
        return PersoGroupMatrix(LGen2, n);
    end;
    INDEF_FORM_AutomorphismGroup:=function(Qmat)
        local n, NSP, TheCompl, FullBasis, p, QmatRed, GRPred, GRPfull, ListGenTot, eGen, eGenB;
        n:=Length(Qmat);
        NSP:=NullspaceIntMat(Qmat);
        TheCompl:=SubspaceCompletionInt(NSP, n);
        FullBasis:=Concatenation(TheCompl, NSP);
        p:=Length(TheCompl);
        QmatRed:=TheCompl * Qmat * TransposedMat(TheCompl);
        GRPred:=INDEF_FORM_AutomorphismGroup_Reduction(QmatRed);
        GRPfull:=ExtendIsometryGroup(GRPred, p, n);
        ListGenTot:=[];
        for eGen in GeneratorsOfGroup(GRPfull)
        do
            eGenB:=Inverse(FullBasis) * eGen * FullBasis;
            if eGenB * Qmat * TransposedMat(eGenB) <> Qmat then
                Error("eGenB should preserve Qmat");
            fi;
            Add(ListGenTot, eGenB);
        od;
        return PersoGroupMatrix(ListGenTot, n);
    end;
    #
    # Vectors in the lattice
    #
    INDEF_FORM_GetOrbitRepresentative:=function(Qmat, X)
        local ApproxModel, ListRepr, FuncInsert, ListCand, eCand;
        Print("Before INDEF_FORM_GetApproximateModel, case 5\n");
        if INDEF_FORM_GetAttackScheme(Qmat).h = 0 then
            return INDEF_FORM_GetOrbitRepresentative_PosNeg(Qmat, X);
        fi;
#        Print("|Qmat|=", Length(Qmat), " X=", X, "\n");
        ApproxModel:=INDEF_FORM_GetApproximateModel(Qmat);
        ListRepr:=[];
        FuncInsert:=function(fRepr)
            local eRepr, test, fNorm;
            fNorm:=fRepr * Qmat * fRepr;
            if fNorm<>X then
                Print("eNorm=", fNorm, " X=", X, "\n");
                Error("The norm is inconsistent");
            fi;
            for eRepr in ListRepr
            do
                test:=INDEF_FORM_EquivalenceVector(Qmat, Qmat, eRepr, fRepr);
                if test<>fail then
                    return;
                fi;
            od;
            Add(ListRepr, fRepr);
        end;
        ListCand:=ApproxModel.GetCoveringOrbitRepresentatives(X);
        for eCand in ListCand
        do
            FuncInsert(eCand);
        od;
        return ListRepr;
    end;
    INDEF_FORM_StabilizerVector:=function(Qmat, v)
        local n, eRec, GRP1, GRP2;
        Print("Beginning of INDEF_FORM_StabilizerVector\n");
        n:=Length(Qmat);
        if RankMat(Qmat) <> n then
            Error("Right now INDEF_FORM_StabilizerVector requires Qmat to be full dimensional");
        fi;
        eRec:=INDEF_FORM_GetVectorStructure(Qmat, v);
        GRP1:=INDEF_FORM_AutomorphismGroup(eRec.GramMatRed);
        GRP2:=eRec.MapOrthogonalSublatticeGroup(GRP1);
        return MatrixIntegral_Stabilizer(n, GRP2);
    end;
    Kernel_Equivalence_Qmat:=function(Qmat1, Qmat2, EquivRat, eRec1, fTest, f_stab)
        local n, GRP1_A, GRP1_B, HasNonIntegralMatrix, eGen, TheRet, ListGen;
        Print("Beginning of Kernel_Equivalence_Qmat\n");
        n:=Length(Qmat1);
        if EquivRat * Qmat1 * TransposedMat(EquivRat) <> Qmat2 then
            Error("EquivRat not mapping Qmat1 to Qmat2");
        fi;
        if IsIntegralMat(EquivRat) then # This is an Ansatz. If the matrix is integral, no need for all the stuff below.
            return EquivRat;
        fi;
        GRP1_A:=f_stab(eRec1);
        Print("We have GRP1_A\n");
        GRP1_B:=eRec1.MapOrthogonalSublatticeGroup(GRP1_A);
        Print("Kernel_Equivalence_Qmat before MatrixIntegral_Equivalence_Bis\n");
        # Find a g1 in GRP1_B such that EquivRat * g1 in GL(n,Z)
        TheRet:=MatrixIntegral_Equivalence_Bis(GRP1_B, EquivRat);
#        Print("Invariant(GRP1_B)=", InvariantSpaceOfGroup(n, GRP1_B), "\n");
        if TheRet = fail then
            return fail;
        fi;
        if TheRet * Qmat1 * TransposedMat(TheRet)<>Qmat2 then
            Error("The redution matrix did not work as we expected it. Please debug");
        fi;
        if fTest(TheRet) = false then
            Error("The fTest has failed");
        fi;
        return TheRet;
    end;
    INDEF_FORM_EquivalenceVector:=function(Qmat1, Qmat2, v1, v2)
        local n, eNorm, eRec1, eRec2, test, EquivRat, GRP1_A, GRP1_B, eGen, TheRet, Subspace1, Subspace2, fTest;
        if INDEF_FORM_InvariantVector(Qmat1, v1) <> INDEF_FORM_InvariantVector(Qmat2, v2) then
            Print("INDEF_FORM_EquivalenceVector returns fail due to different invariants\n");
            return fail;
        fi;
        Print("INDEF_FORM_EquivalenceVector v1=", v1, " v2=", v2, "\n");
        n:=Length(Qmat1);
        if RankMat(Qmat1)<>n then
            Error("Right now INDEF_FORM_EquivalenceVector requires Qmat1 to be full dimensional");
        fi;
        eNorm:=v1 * Qmat1 * v1;
        eRec1:=INDEF_FORM_GetVectorStructure(Qmat1, v1);
        eRec2:=INDEF_FORM_GetVectorStructure(Qmat2, v2);
        test:=INDEF_FORM_TestEquivalence(eRec1.GramMatRed, eRec2.GramMatRed);
        if test=fail then
            Print("INDEF_FORM_EquivalenceVector returns fail because of not equivalent GramMatRed\n");
            return fail;
        fi;
        if eNorm<>0 then
            EquivRat:=eRec2.PmatInv * ExpandMatrix(test) * eRec1.Pmat;
        else
            Subspace1:=Inverse(test) * eRec2.NSP;
            Subspace2:=eRec1.NSP;
            Print("Before calling LORENTZ_ExtendOrthogonalIsotropicIsomorphism, case 2\n");
            EquivRat:=LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(Qmat1, Subspace1, Qmat2, Subspace2);
            if v1 * Inverse(EquivRat) = -v2 then
                EquivRat:=-EquivRat;
            fi;
        fi;
        if v1*Inverse(EquivRat)<>v2 then
            Error("Vectors are not mapped correctly");
        fi;
        fTest:=function(eTest)
            return v1*Inverse(eTest) = v2;
        end;
        return Kernel_Equivalence_Qmat(Qmat1, Qmat2, EquivRat, eRec1, fTest, f_stab_plane);
    end;
    #
    # Isotropic k-planes
    # stuff = "plane" or "flag"
    f_equiv_plane:=function(eRec1, eRec2)
        return INDEF_FORM_TestEquivalence(eRec1.GramMatRed, eRec2.GramMatRed);
    end;
    f_stab_plane:=function(eRec)
        return INDEF_FORM_AutomorphismGroup(eRec.GramMatRed);
    end;
    INDEF_FORM_Equivalence_IsotropicKstuff_Kernel:=function(Qmat1, Qmat2, Plane1, Plane2, f_equiv, f_stab)
        local eRec1, eRec2, test, Subspace1, Subspace2, EquivRat, fTest, TheRec;
        Print("Beginning of INDEF_FORM_Equivalence_IsotropicKplane\n");
        if INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat1, Plane1, f_stab) <> INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat2, Plane2, f_stab) then
            Print("INDEF_FORM_Equivalence_IsotropicKplane returns fail due to different invariants\n");
            return fail;
        fi;
        eRec1:=INDEF_FORM_GetRec_IsotropicKplane(Qmat1, Plane1);
        eRec2:=INDEF_FORM_GetRec_IsotropicKplane(Qmat2, Plane2);
        test:=f_equiv(eRec1, eRec2);
        if test=fail then
            Print("INDEF_FORM_EquivalenceVector returns fail because of not equivalent GramMatRed\n");
            return fail;
        fi;
        Subspace1:=Inverse(test) * eRec2.NSP;
        Subspace2:=eRec1.NSP;
        Print("Before calling LORENTZ_ExtendOrthogonalIsotropicIsomorphism, case 2\n");
        TheRec:=LORENTZ_ExtendOrthogonalIsotropicIsomorphism(Qmat1, Subspace1, Qmat2, Subspace2);
        EquivRat:=TheRec.get_one_transformation();
        if TestEqualitySpace(Plane1 * Inverse(EquivRat), Plane2) = false then
            Error("Plane1 and Plane2 should be mapped");
        fi;
        fTest:=function(eTest)
            return TestEqualitySpace(Plane1*Inverse(eTest), Plane2);
        end;
        return Kernel_Equivalence_Qmat(Qmat1, Qmat2, EquivRat, eRec1, fTest, f_stab);
    end;
    INDEF_FORM_Equivalence_IsotropicKplane:=function(Qmat1, Qmat2, Plane1, Plane2)
        return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(Qmat1, Qmat2, Plane1, Plane2, f_equiv_plane, f_stab_plane);
    end;
    INDEF_FORM_Invariant_IsotropicKstuff_Kernel:=function(Qmat, Plane, f_stab)
        local n, eRec, GRP1, GRP2, eInvRed;
        n:=Length(Qmat);
        eRec:=INDEF_FORM_GetRec_IsotropicKplane(Qmat, Plane);
        GRP1:=f_stab(eRec);
        GRP2:=eRec.MapOrthogonalSublatticeGroup(GRP1);
        eInvRed:=INDEF_FORM_Invariant_IsotropicKplane_Raw(Qmat, Plane);
        return [GetRationalInvariant(GRP2), eInvRed];
    end;
    INDEF_FORM_Invariant_IsotropicKplane:=function(Qmat, Plane)
        return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat, Plane, f_stab_plane);
    end;
    INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel:=function(Qmat, Plane, f_stab)
        local n, eRec, GRP1, GRP2;
        Print("Beginning of INDEF_FORM_StabilizerVector\n");
        n:=Length(Qmat);
        if RankMat(Qmat) <> n then
            Error("Right now INDEF_FORM_StabilizerVector requires Qmat to be full dimensional");
        fi;
        eRec:=INDEF_FORM_GetRec_IsotropicKplane(Qmat, Plane);
        GRP1:=f_stab(eRec);
        GRP2:=eRec.MapOrthogonalSublatticeGroup(GRP1);
        return MatrixIntegral_Stabilizer(n, GRP2);
    end;
    INDEF_FORM_Stabilizer_IsotropicKplane:=function(Qmat, Plane)
        return INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(Qmat, Plane, f_stab_plane);
    end;
    INDEF_FORM_RightCosets_IsotropicKstuff_Kernel:=function(Qmat, ePlane, f_stab)
        local n, eRec, GRP1, GRP2, ListRightCoset, eGen, ePlaneImg, eLine, eCos;
        # We have two groups:
        # -- The group stabilizing ePlane
        # -- The group stabilizing ePlane^{perp} and its mapping to the full group.
        # The group stabilizing ePlane^{perp} can be injected 
        Print("Beginning of INDEF_FORM_StabilizerVector\n");
        n:=Length(Qmat);
        if RankMat(Qmat) <> n then
            Error("Right now INDEF_FORM_StabilizerVector requires Qmat to be full dimensional");
        fi;
        eRec:=INDEF_FORM_GetRec_IsotropicKplane(Qmat, ePlane);
        GRP1:=f_stab(eRec);
        for eGen in GeneratorsOfGroup(GRP1)
        do
            if eGen * eRec.GramMatRed * TransposedMat(eGen) <> eRec.GramMatRed then
                Error("GRP1 does not preserve eRec.GramMatRed");
            fi;
        od;
        GRP2:=eRec.MapOrthogonalSublatticeGroup(GRP1);
        for eGen in GeneratorsOfGroup(GRP2)
        do
            if eGen * Qmat * TransposedMat(eGen) <> Qmat then
                Error("GRP2 does not preserve Qmat");
            fi;
            ePlaneImg:=ePlane * eGen;
            for eLine in ePlaneImg
            do
                if SolutionIntMat(ePlane, eLine)=fail then
                    Error("ePlane should be left invariant by eGen");
                fi;
            od;
        od;
        ListRightCoset:=MatrixIntegral_RightCosets(n, GRP2);
        for eCos in ListRightCoset
        do
            if eCos * Qmat * TransposedMat(eCos)<>Qmat then
                Error("eCos is not preserving Qmat");
            fi;
        od;
        return ListRightCoset;
    end;
    INDEF_FORM_RightCosets_IsotropicKplane:=function(Qmat, ePlane)
        return INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(Qmat, ePlane, f_stab_plane);
    end;
    INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel:=function(Qmat, k, RecF)
        local eNorm, ListOrbit, iK, ListRecReprKplane, fInsert, fGetRepresentatives, eRepr, iRepr, eRecReprExp;
        eNorm:=0;
        ListOrbit:=List(INDEF_FORM_GetOrbitRepresentative(Qmat, eNorm), x->[x]);
        Print("k=1 ListOrbit=", ListOrbit, "\n");
        for iK in [2..k]
        do
            ListRecReprKplane:=[];
            fInsert:=function(fRecReprKplane)
                local eRecReprKplane, test;
                for eRecReprKplane in ListRecReprKplane
                do
                    if eRecReprKplane.eInv=fRecReprKplane.eInv then
                        test:=RecF.f_equiv(Qmat, Qmat, eRecReprKplane.ePlane, fRecReprKplane.ePlane);
                        if test<>fail then
                            return;
                        fi;
                    fi;
                od;
                Add(ListRecReprKplane, fRecReprKplane);
                Print("iK=", iK, " now |ListRecReprKplane|=", Length(ListRecReprKplane), "\n");
            end;
            fGetRepresentatives:=function(Qmat, ePlane)
                local NSP, dimNSP, ePlane_expr, eV, eSol, ComplBasisInNSP, NSP_sub, QmatRed, ListOrbitF, ListRecReprRet, eVect, eVectB, eVectC, ListRightCosets, eCos, ePlaneB, eInv;
                # Some possible improvement. Use the double cosets
                # The double coset consists in splitting an orbit x G as
                # y1 H \cup ..... yN H
                # Or in other words G = \cup_i Stab(x) y_i H
                #
                # In that case
                # G = group of rational transformation preserving L = S^{perp}
                #     and acting integrally on it.
                # Stab(x) = Integral transformations preserving L and some vector x in L.
                # H = group of integral transformations preserving the big lattice
                #     and preserving L.
                #
                # If G = H the only one entry to treat.
                # What that mean is that when we go down the chain of subgroup by breaking down
                # the orbit split. Can this happen in the same way over all the orbits?
                # The full of the lattice group is G(Qmat) and the full orbit is x G(Qmat)
                # We want to write the code as 
                NSP:=NullspaceIntMat(TransposedMat(ePlane * Qmat));
                dimNSP:=Length(NSP);
                ePlane_expr:=[];
                for eV in ePlane
                do
                    eSol:=SolutionIntMat(NSP, eV);
                    if eSol=fail then
                        Error("eSol should not be fail");
                    fi;
                    Add(ePlane_expr, eSol);
                od;
                ComplBasisInNSP:=SubspaceCompletionInt(ePlane_expr, dimNSP);
                NSP_sub:=ComplBasisInNSP * NSP;
                QmatRed:=NSP_sub * Qmat * TransposedMat(NSP_sub);
                ListOrbitF:=INDEF_FORM_GetOrbitRepresentative(QmatRed, eNorm);
                Print("|ListOrbitF|=", Length(ListOrbitF), "\n");
                ListRightCosets:=RecF.f_coset(Qmat, ePlane);
                Print("|ListRightCosets|=", Length(ListRightCosets), "\n");
                ListRecReprRet:=[];
                for eVect in ListOrbitF
                do
                    eVectB:=eVect * NSP_sub;
                    for eCos in ListRightCosets
                    do
                        eVectC:=eVectB * eCos;
                        if eVectC * Qmat * eVectC <> 0 then
                            Error("eVectC is not isotropic");
                        fi;
                        ePlaneB:=Concatenation(ePlane, [eVectC]);
                        if IsIntegralMat(ePlaneB)=false then
                            Error("The matrix should be integral");
                        fi;
                        eInv:=RecF.f_inv(Qmat, ePlaneB);
                        Add(ListRecReprRet, rec(eInv:=eInv, ePlane:=ePlaneB));
                    od;
                od;
                return ListRecReprRet;
            end;
            iRepr:=0;
            for eRepr in ListOrbit
            do
                iRepr:=iRepr+1;
                Print("At iK=", iK, " iRepr=", iRepr, " eRepr=", eRepr, "\n");
                for eRecReprExp in fGetRepresentatives(Qmat, eRepr)
                do
                    fInsert(eRecReprExp);
                od;
            od;
            ListOrbit:=List(ListRecReprKplane, x->x.ePlane);
        od;
        return ListOrbit;
    end;
    INDEF_FORM_GetOrbit_IsotropicKplane:=function(Qmat, k)
        local RecF;
        RecF:=rec(f_equiv:=INDEF_FORM_Equivalence_IsotropicKplane,
                  f_coset:=INDEF_FORM_RightCosets_IsotropicKplane,
                  f_inv:=INDEF_FORM_Invariant_IsotropicKplane);
        return INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(Qmat, k, RecF);
    end;
    #
    # Isotropic k-flags
    #
    f_equiv_flag:=function(eRec1, eRec2)
        local PlaneExpr1, PlaneExpr2, the_dim, TheCompl1, TheCompl2, FullBasis1, FullBasis2, QmatRed1, QmatRed2, test, TheEquivTest, p, i, j, TheEquiv, ListSpaces1, ListSpaces2;
        PlaneExpr1:=List(eRec1.Plane, x->SolutionIntMat(eRec1.NSP, x));
        PlaneExpr2:=List(eRec2.Plane, x->SolutionIntMat(eRec2.NSP, x));
        the_dim:=Length(eRec1.NSP);
        TheCompl1:=SubspaceCompletionInt(PlaneExpr1, the_dim);
        TheCompl2:=SubspaceCompletionInt(PlaneExpr2, the_dim);
        FullBasis1:=Concatenation(TheCompl1, PlaneExpr1);
        FullBasis2:=Concatenation(TheCompl2, PlaneExpr2);
        QmatRed1:=TheCompl1 * eRec1.GramMatRed * TransposedMat(TheCompl1);
        QmatRed2:=TheCompl2 * eRec2.GramMatRed * TransposedMat(TheCompl2);
        test:=INDEF_FORM_TestEquivalence(QmatRed1, QmatRed2);
        if test=fail then
            return fail;
        fi;
        TheEquivTest:=IdentityMat(the_dim);
        p:=Length(QmatRed1);
        for i in [1..p]
        do
            for j in [1..p]
            do
                TheEquivTest[i][j]:=test[i][j];
            od;
        od;
        TheEquiv:=Inverse(FullBasis2) * TheEquivTest * FullBasis1;
        if TheEquiv * eRec1.GramMatRed * TransposedMat(TheEquiv) <> eRec2.GramMatRed then
            Error("This is not an equivalence");
        fi;
        ListSpaces1:=f_get_list_spaces(PlaneExpr1);
        ListSpaces2:=f_get_list_spaces(PlaneExpr2);
        for i in [1..Length(ListSpaces1)]
        do
            if TestEqualitySpace(ListSpaces1[i] * Inverse(TheEquiv), ListSpaces2[i]) = false then
                Error("The space are not correctly mapped");
            fi;
        od;
        return TheEquiv;
    end;
    f_stab_flag:=function(eRec)
        local k, PlaneExpr, the_dim, TheCompl, FullBasis, QmatRed, GRPred, GRPfull, ListSpaces, ListGenTot, eGen, eGenB, eSpace;
        k:=Length(eRec.Plane);
        PlaneExpr:=List(eRec.Plane, x->SolutionIntMat(eRec.NSP, x));
        the_dim:=Length(eRec.NSP);
        TheCompl:=SubspaceCompletionInt(PlaneExpr, the_dim);
        FullBasis:=Concatenation(TheCompl, PlaneExpr);
        QmatRed:=TheCompl * eRec.GramMatRed * TransposedMat(TheCompl);
        GRPred:=INDEF_FORM_AutomorphismGroup(QmatRed);
        GRPfull:=ExtendIsometryGroup_Triangular(GRPred, Length(TheCompl), the_dim);
        ListSpaces:=f_get_list_spaces(PlaneExpr);
        ListGenTot:=[];
        for eGen in GeneratorsOfGroup(GRPfull)
        do
            eGenB:=Inverse(FullBasis) * eGen * FullBasis;
            if eGenB * eRec.GramMatRed * TransposedMat(eGenB) <> eRec.GramMatRed then
                Error("eGenB should preserve eRec.GramMatRed");
            fi;
            for eSpace in ListSpaces
            do
                if TestEqualitySpace(eSpace * eGenB, eSpace) = false then
                    Error("The space eSpace is not correctly preserved");
                fi;
            od;
            Add(ListGenTot, eGenB);
        od;
        return PersoGroupMatrix(ListGenTot, the_dim);
    end;
    INDEF_FORM_Equivalence_IsotropicKflag:=function(Qmat1, Qmat2, Plane1, Plane2)
        return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(Qmat1, Qmat2, Plane1, Plane2, f_equiv_flag, f_stab_flag);
    end;
    INDEF_FORM_Invariant_IsotropicKflag:=function(Qmat, Plane)
        return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat, Plane, f_stab_flag);
    end;
    INDEF_FORM_Stabilizer_IsotropicKflag:=function(Qmat, Plane)
        return INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(Qmat, Plane, f_stab_flag);
    end;
    INDEF_FORM_RightCosets_IsotropicKflag:=function(Qmat, ePlane)
        return INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(Qmat, ePlane, f_stab_flag);
    end;
    INDEF_FORM_GetOrbit_IsotropicKflag:=function(Qmat, k)
        local RecF;
        RecF:=rec(f_equiv:=INDEF_FORM_Equivalence_IsotropicKflag,
                  f_coset:=INDEF_FORM_RightCosets_IsotropicKflag,
                  f_inv:=INDEF_FORM_Invariant_IsotropicKflag);
        return INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(Qmat, k, RecF);
    end;
    return rec(INDEF_FORM_GetOrbitRepresentative:=INDEF_FORM_GetOrbitRepresentative,
               INDEF_FORM_StabilizerVector:=INDEF_FORM_StabilizerVector,
               INDEF_FORM_EquivalenceVector:=INDEF_FORM_EquivalenceVector,
               INDEF_FORM_StabilizerVector:=INDEF_FORM_StabilizerVector,
               INDEF_FORM_EquivalenceVector:=INDEF_FORM_EquivalenceVector,
               INDEF_FORM_AutomorphismGroup:=INDEF_FORM_AutomorphismGroup,
               INDEF_FORM_TestEquivalence:=INDEF_FORM_TestEquivalence,
               INDEF_FORM_Invariant_IsotropicKplane:=INDEF_FORM_Invariant_IsotropicKplane,
               INDEF_FORM_Equivalence_IsotropicKplane:=INDEF_FORM_Equivalence_IsotropicKplane,
               INDEF_FORM_Stabilizer_IsotropicKplane:=INDEF_FORM_Stabilizer_IsotropicKplane,
               INDEF_FORM_RightCosets_IsotropicKplane:=INDEF_FORM_RightCosets_IsotropicKplane,
               INDEF_FORM_GetOrbit_IsotropicKplane:=INDEF_FORM_GetOrbit_IsotropicKplane,
               INDEF_FORM_Equivalence_IsotropicKflag:=INDEF_FORM_Equivalence_IsotropicKflag,
               INDEF_FORM_Invariant_IsotropicKflag:=INDEF_FORM_Invariant_IsotropicKflag,
               INDEF_FORM_Stabilizer_IsotropicKflag:=INDEF_FORM_Stabilizer_IsotropicKflag,
               INDEF_FORM_RightCosets_IsotropicKflag:=INDEF_FORM_RightCosets_IsotropicKflag,
               INDEF_FORM_GetOrbit_IsotropicKflag:=INDEF_FORM_GetOrbit_IsotropicKflag);
end;


INDEF_FORM_GetOrbitRepresentative:=function(Qmat, X_v)
    local X;
    X:=X_v[1];
    return INDEF_FORM_Machinery_AllFct().INDEF_FORM_GetOrbitRepresentative(Qmat,X);
end;


INDEF_FORM_StabilizerVector:=function(Qmat, v)
    return INDEF_FORM_Machinery_AllFct().INDEF_FORM_StabilizerVector(Qmat, v);
end;


INDEF_FORM_EquivalenceVector:=function(Qmat1, Qmat2, v1, v2)
    return INDEF_FORM_Machinery_AllFct().INDEF_FORM_EquivalenceVector(Qmat1, Qmat2, v1, v2);
end;


INDEF_FORM_AutomorphismGroup:=function(Qmat)
    return INDEF_FORM_Machinery_AllFct().INDEF_FORM_AutomorphismGroup(Qmat);
end;


INDEF_FORM_TestEquivalence:=function(Qmat1, Qmat2)
    return INDEF_FORM_Machinery_AllFct().INDEF_FORM_TestEquivalence(Qmat1, Qmat2);
end;
