CVPVallentinProgram_Rational:=function(GramMat, eV)
    local GramMat_oscar, eV_oscar, ListVect_oscar;
    GramMat_oscar:=MatrixToOscar(GramMat);
    eV_oscar:=VectorToOscar(eV);
    ListVect_oscar:=Julia.Indefinite.LATT_near(GramMat_oscar, eV_oscar, 0);
    return ReadOscarMatrix(ListVect_oscar);
end;

ClosestAtDistanceVallentinProgram:=function(GramMat, eV, TheDist)
    local GramMat_oscar, eV_oscar, TheDist_oscar, ListVect_oscar;
    GramMat_oscar:=MatrixToOscar(GramMat);
    eV_oscar:=VectorToOscar(eV);
    TheDist_oscar:=ScalarToOscar(TheDist);
    ListVect_oscar:=Julia.Indefinite.LATT_near(GramMat_oscar, eV_oscar, TheDist_oscar);
    return ReadOscarMatrix(ListVect_oscar);
end;
