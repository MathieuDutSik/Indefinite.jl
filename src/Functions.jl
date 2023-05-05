function INDEF_FORM_GetOrbitRepresentative(Qmat::Nemo.QQMatrix, Xval::Nemo.QQFieldElem)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  Xval_gap = GAP.evalstr(string("[", string(Xval), "]"))
  print("typeof(Qmat_gap)=", typeof(Qmat_gap), "\n")
  print("typeof(Xval_gap)=", typeof(Xval_gap), "\n")
  ListOrbitRepr_gap = GAP.Globals.INDEF_FORM_GetOrbitRepresentative(Qmat_gap, Xval_gap)
  return GAP.Globals.MatrixToOscar(ListOrbitRepr_gap)
end

function INDEF_FORM_AutomorphismGroup(Qmat::Nemo.QQMatrix)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  GRP_gap = GAP.Globals.INDEF_FORM_AutomorphismGroup(Qmat_gap)
  ListGens_gap = GAP.Globals.GeneratorsOfGroup(GRP_gap)
  return GAP.Globals.ListMatrixToOscar(ListGens_gap)
end