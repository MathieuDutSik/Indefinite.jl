ScalarToOscarString:=function(eScal)
    local eNum, eDen, eStr;
    if IsInt(eScal) then
        return String(eScal);
    else
        eNum:=NumeratorRat(eScal);
        eDen:=DenominatorRat(eScal);
        eStr:=Concatenation(String(eNum), "//", String(eDen));
        return eStr;
    fi;
end;


MatrixToOscarString:=function(M)
    local nbLine, nbCol, TheStr, IsFirst, eLine, eVal;
    nbLine:=Length(M);
    if nbLine = 0 then
        return "matrix(Nemo.QQ, 0, 0)";
    fi;
    nbCol:=Length(M[1]);
    TheStr:=Concatenation("matrix(Nemo.QQ,", String(nbLine), ",", String(nbCol), ",[");
    IsFirst:=true;
    for eLine in M
    do
        for eVal in eLine
        do
            if IsFirst=false then
                TheStr:=Concatenation(TheStr, ",");
            fi;
            IsFirst:=false;
            TheStr:=Concatenation(TheStr, ScalarToOscarString(eVal));
        od;
    od;
    TheStr:=Concatenation(TheStr, "])");
    return TheStr;
end;


MatrixToOscar:=function(M)
    local M_str;
    M_str:=MatrixToOscarString(M);
    Print("M_str=", M_str, "\n");
    return JuliaEvalString(M_str);
end;


ListMatrixToOscar:=function(ListM)
    local TotStr, i;
    Print("LustM=", ListM, "\n");
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


VectorToOscar:=function(V)
    return MatrixToOscar([V]);
end;


StringToOscar:=function(estr)
    local fstr;
    fstr:=Concatenation("\"", estr, "\"");
    return JuliaEvalString(fstr);
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
    local M;
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
