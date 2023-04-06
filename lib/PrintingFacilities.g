

PARI_PrintMatrix:=function(output, eMat)
  local nbLine, nbCol, iLine, iCol;
  nbLine:=Length(eMat);
  nbCol:=Length(eMat[1]);
  if nbLine=1 and nbCol=1 then
    AppendTo(output, "matrix(1,1,i,j,1)");
    return;
  fi;
  AppendTo(output, "[");
  for iLine in [1..nbLine]
  do
    if iLine>1 then
      AppendTo(output, ";");
    fi;
    for iCol in [1..nbCol]
    do
      if iCol>1 then
        AppendTo(output, ",");
      fi;
      AppendTo(output, eMat[iLine][iCol]);
    od;
  od;
  AppendTo(output, "]");
end;
