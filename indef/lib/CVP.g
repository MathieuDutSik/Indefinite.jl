ClosestAtDistanceVallentinProgram:=function(GramMat, eV, TheDist)
    local GramMat_oscar, eV_oscar, TheDist_oscar, ListVect_oscar;
    GramMat_oscar:=MatrixToOscar(GramMat);
    eV_oscar:=VectorToOscar(eV);
    TheDist_oscar:=MatrixToOscar([[TheDist]]);
    ListVect_oscar:=Julia.Indefinite.LATT_near(GramMat_oscar, eV_oscar, TheDist_oscar);
    return ReadOscarMatrix(ListVect_oscar);
end;

CVPVallentinProgram:=function(GramMat, eV)
    return ClosestAtDistanceVallentinProgram(GramMat, eV, 0);
end;

