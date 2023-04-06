

ShortVectorDutourVersion:=function(GramMat, eNorm)
  local H, eNormRed, GramMatRed, SHV, ListVector, iV, eV;
  H:=RemoveFractionMatrixPlusCoef(GramMat);
  eNormRed:=LowerInteger(eNorm*H.TheMult);
  GramMatRed:=H.TheMat;
  if IsPositiveDefiniteSymmetricMatrix(GramMatRed)=false then
    Error("Matrix should be positive definite");
  fi;
  SHV:=ShortestVectors(GramMatRed, eNormRed);
  ListVector:=[];
  for iV in [1..Length(SHV.vectors)]
  do
    if SHV.norms[iV] <= eNormRed then
      eV:=SHV.vectors[iV];
      Add(ListVector, eV);
      Add(ListVector, -eV);
    fi;
  od;
  return ListVector;
end;
