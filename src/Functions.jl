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

function INDEF_FORM_TestEquivalence(Qmat1::Nemo.QQMatrix, Qmat2::Nemo.QQMatrix)
  Qmat1_gap = GAP.Globals.ReadOscarMatrix(Qmat1)
  Qmat2_gap = GAP.Globals.ReadOscarMatrix(Qmat2)
  TheEquiv_gap = GAP.Globals.INDEF_FORM_TestEquivalence(Qmat1_gap, Qmat2_gap)
  return GAP.Globals.MatrixToOscar(TheEquiv_gap)
end

function INDEF_FORM_GetOrbit_IsotropicKplane(Qmat::Nemo.QQMatrix, k::Int64)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  ListRepr_gap = GAP.Globals.INDEF_FORM_GetOrbit_IsotropicKplane(Qmat_gap, k)
  return GAP.Globals.ListMatrixToOscar(ListRepr_gap)
end

function INDEF_FORM_GetOrbit_IsotropicKflag(Qmat::Nemo.QQMatrix, k::Int64)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  ListRepr_gap = GAP.Globals.INDEF_FORM_GetOrbit_IsotropicKflag(Qmat_gap, k)
  return GAP.Globals.ListMatrixToOscar(ListRepr_gap)
end

