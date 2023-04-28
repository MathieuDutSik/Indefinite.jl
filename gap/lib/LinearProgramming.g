#
# this thing is the end of a long design road.
# realizing that linear programming is quite complex.
# Basically, everything is outputed, everything is read
# and you have to make interpretations yourself.
LinearProgramming:=function(InequalitySet, ToBeMinimized)
    InequalitySet_oscar:=MatrixToOscar(InequalitySet);
    ToBeMinimized_oscar:=VectorToOscar(ToBeMinimized);
    #
    TheResult:=Oscar.LinearProgramming(InequalitySet_oscar, ToBeMinimized_oscar);
    answer:=TheResult[1];
    DirectSolution:=OscarVectorToVector(TheResult[2]);
    DualSolution:=OscarVectorToVector(TheResult[3]);
    return rec(DirectSolution:=DirectSolution, DualSolution:=DualSolution);
end;


GetInitialRaysGeneral:=function(FAC, command)
    FAC_oscar:=MatrixToOscar(FAC);
    TheResult:=Oscar.POLY_samplingFacets(FAC_oscar, command);
    return ReadOscarListIncd(TheResult);
end;

GetInitialRays_LinProg:=function(EXT, nb)
    local command;
    command:=Concatenation("lp_cdd:iter_", String(nb));
    return GetInitialRaysGeneral(EXT, command);
end;


Norm_L1:=function(eVect)
  return Sum(List(eVect, AbsInt));
end;
