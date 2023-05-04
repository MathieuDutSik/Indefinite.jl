


function INDEF_FORM_GetOrbitRepresentative(Qmat::Nemo.QQMatrix, Xval::Nemo.QQFieldElem)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  Xval_gap = GAP.evalstr(string(Xval))
  ListOrbitRepr_gap = INDEF_FORM_GetOrbitRepresentative(Qmat_gap, Xval_gap)
  return GAP.Globals.MatrixToOscar(ListOrbitRepr_gap)
end
