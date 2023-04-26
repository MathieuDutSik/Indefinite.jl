MatrixToOscarString:=function(M)
    local TheStr, IsFirst, eLine, eVal;
    TheStr:="[";
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
    TheStr:=Concatenation(TheStr, "]");
    return TheStr;
end;

MatrixToOscar:=function(M)
    local M_str;
    M_str:=MatrixToOscarString(M);
    return JuliaEvalString(M_str);
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
    local M_str;
    V_str:=VectorToOscarString(V);
    return JuliaEvalString(V_str);
end;

ReadOscarVector:=function(V_oscar)
    local eList;
    eList:=JuliaToGAP(IsList, V_oscar);
    return List(eList, Oscar.GAP.julia_to_gap);
end;

ReadOscarListIncd:=function(ListIncd_oscar)
    eList:=JuliaToGAP(IsList, V_oscar);
    TheList:=[];
    for eEnt in eList
    do
        Add(TheList, OscarVectorToVector(eEnt));
    od;
    return TheList;
end;


PermutationGroupToOscar:=function(n, PermGroup)
    local ListGen, TheStr, IsFirst, eGen, i, eImg, n_gen, PermGRoup_oscar;
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
