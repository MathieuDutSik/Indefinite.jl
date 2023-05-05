#
# this thing is the end of a long design road.
# realizing that linear programming is quite complex.
# Basically, everything is outputed, everything is read
# and you have to make interpretations yourself.
LinearProgramming:=function(InequalitySet, ToBeMinimized)
    local InequalitySet_oscar, ToBeMinimized_oscar, TheResult, answer, optimal_value, DirectSolution, DualSolution;
    InequalitySet_oscar:=MatrixToOscar(InequalitySet);
    ToBeMinimized_oscar:=VectorToOscar(ToBeMinimized);
    #
    TheResult:=Julia.Indefinite.LinearProgramming(InequalitySet_oscar, ToBeMinimized_oscar);
    answer:=TheResult[1];
    optimal_value:=TheResult[2];
    DirectSolution:=ReadOscarVector(TheResult[2]);
    DualSolution:=ReadOscarVector(TheResult[3]);
    return rec(answer:=answer, optimalValue:=optimal_value, DirectSolution:=DirectSolution, DualSolution:=DualSolution);
end;


GetInitialRaysGeneral:=function(FAC, command)
    local FAC_oscar, TheResult, command_oscar;
    FAC_oscar:=MatrixToOscar(FAC);
    command_oscar:=StringToOscar(command);
    TheResult:=Julia.Indefinite.POLY_sampling_facets(FAC_oscar, command_oscar);
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
