ArithmeticAutomorphismGroup:=function(ListGramMat)
    local ListGramMat_oscar, ListGens_oscar, ListGens, eGen, eGram;
    ListGramMat_oscar:=ListMatrixToOscar(ListGramMat);
    ListGens_oscar:=Oscar.LATT_Automorphism(ListGramMat_oscar);
    ListGens:=ReadOscarListMatrix(ListGens_oscar);
    for eGen in ListGens
    do
        for eGram in ListGramMat
        do
            if eGram<>eGen*eGram*TransposedMat(eGen) then
                Error("The Gram matrix should be preserved");
            fi;
        od;
    od;
    return Group(ListGens);
end;


ArithmeticIsomorphism:=function(ListGramMat1, ListGramMat2)
    local ListGramMat1_oscar, ListGramMat2_oscar, TheMat_oscar, TheMat, i, eGram1, eGram2;
    ListGramMat1_oscar:=ListMatrixToOscar(ListGramMat1);
    ListGramMat2_oscar:=ListMatrixToOscar(ListGramMat2);
    TheMat_oscar:=Oscar.LATT_Isomorphism(ListGramMat1_oscar, ListGramMat2_oscar);
    TheMat:=ReadOscarMatrix(TheMat_oscar);
    if Length(TheMat)=0 then
        return false;
    fi;
    for i in [1..Length(ListGramMat1)]
    do
        eGram1:=ListGramMat1[i];
        eGram2:=ListGramMat2[i];
        if TheMat * eGram1 * TransposedMat(TheMat) <> eGram2 then
            Error("The eGram1 is not mapped to eGram2");
        fi;
    od;
    return TheMat;
end;
