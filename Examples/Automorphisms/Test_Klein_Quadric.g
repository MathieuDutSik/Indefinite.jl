Hmat:=[[0,1],[1,0]];
LorMat:=LORENTZ_AssembleDiagBlocks([Hmat, Hmat, Hmat]);
n:=Length(LorMat);


eFile:="DATA/Klein_Group_33";
if IsExistingFile(eFile) then
    OrthogonalGroup:=ReadAsFunction(eFile)();
else
    OrthogonalGroup:=INDEF_FORM_AutomorphismGroup(LorMat);
    SaveDataToFile(eFile, OrthogonalGroup);
fi;
