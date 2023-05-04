MatrixToOscarString:=function(M)
    local TheStr, IsFirst, eLine, eVal;
    nbLine:=Length(M);
    if nbLine = 0 then
        return "matrix(Nemo.QQ, 0, 0)";
    fi;
    nbCol:=Length(M[1]);
    TheStr:=Concatenation("matrix(QQ,", String(nbLine), ",", String(nbCol), ",[";
    IsFirst:=true;
    for eLine in M
    do
        for eVal in eLine
        do
            if IsFirst=false then
                TheStr:=Concatenation(TheStr, ",");
            fi;
            IsFirst:=false;
            TheStr:=Concatenation(TheStr, String(eVal));
        od;
    od;
    TheStr:=Concatenation(TheStr, "])");
    return TheStr;
end;

MatrixToOscar:=function(M)
    local M_str;
    M_str:=MatrixToOscarString(M);
    return JuliaEvalString(M_str);
end;


ListMatrixToOscar:=function(ListM)
    local TotStr, i;
    TotStr:="Vector{QQMatrix}([";
    for i in [1..Length(ListM)]
    do
        if i > 1 then
            TotStr:=Concatenation(TotStr, ",");
        fi;
        TotStr:=Concatenation(TotStr, MatrixToOscarString(ListM[i]));
    od;
    TotStr:=Concatenation(TotStr, "])");
    return JuliaEvalString(TotStr);
end;



ScalarToOscar:=function(eScal)
    return Oscar.Parse_QQ(String(eScal));
end;


VectorToOscarString:=function(eVect)
    local TheStr, IsFirst, eLine, eVal;
    TheStr:="[";
    IsFirst:=true;
    for eVal in eVect
    do
        if IsFirst=false then
            TheStr:=Concatenation(TheStr, ",");
        fi;
        IsFirst:=false;
        TheStr:=Concatenation(TheStr, String(eVal));
    od;
    TheStr:=Concatenation(TheStr, "]");
    return TheStr;
end;

VectorToOscar:=function(V)
    local V_str;
    V_str:=VectorToOscarString(V);
    return JuliaEvalString(V_str);
end;

ReadOscarMatrix:=function(M_oscar)
    local nbLine, nbCol, M, iLine, eLine, iCol, val;
    Print("M_oscar=", M_oscar, "\n");
    nbLine:=Oscar.Nemo.nrows(M_oscar);
    if nbLine=0 then
        return [];
    fi;
    nbCol:=Oscar.Nemo.ncols(M_oscar);
    M:=[];
    for iLine in [1..nbLine]
    do
        eLine:=[];
        for iCol in [1..nbCol]
        do
            val:=Oscar.GAP.julia_to_gap(M_oscar[iLine,iCol]);
            Add(eLine, val);
        od;
        Add(M, eLine);
    od;
    Print("M=", M, "\n");
    return M;
end;

ReadOscarVector:=function(V_oscar)
    local eList;
    M:=ReadOscarMatrix(V_oscar);
    return M[1];
end;

ReadOscarListIncd:=function(ListIncd_oscar)
    local TheMat, ListIncd, eLine, n, eIncd;
    TheMat:=ReadOscarMatrix(ListIncd_oscar);
    ListIncd:=[];
    for eLine in TheMat
    do
        n:=Length(eLine);
        eIncd:=Filtered([1..n], x->eLine[x]=1);
        Add(ListIncd, eIncd);
    od;
    return ListIncd;
end;

ReadOscarListMatrix:=function(ListMatrix_oscar)
    local eList;
    eList:=JuliaToGAP(IsList, ListMatrix_oscar);
    return List(eList, ReadOscarMatrix);
end;




PermutationGroupToOscar:=function(n, PermGroup)
    local ListGen, TheStr, IsFirst, eGen, i, eImg, n_gen, PermGroup_oscar;
    ListGen:=GeneratorsOfGroup(PermGroup);
    TheStr:="[";
    IsFirst:=true;
    for eGen in ListGen
    do
        for i in [1..n]
        do
            if IsFirst=false then
                TheStr:=Concatenation(TheStr, ",");
            fi;
            IsFirst:=false;
            eImg:=OnPoints(i, eGen);
            TheStr:=Concatenation(TheStr, String(eImg));
        od;
    od;
    TheStr:=Concatenation(TheStr, "]");
    n_gen:=Length(ListGen);
    PermGroup_oscar := Oscar.matrix(Oscar.ZZ, n, n_gen, JuliaEvalString(TheStr));
    return PermGroup_oscar;
end;
