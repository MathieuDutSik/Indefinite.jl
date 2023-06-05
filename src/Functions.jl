"""
   INDEF_FORM_GetOrbitRepresentative(Qmat::Nemo.QQMatrix, Xval::Nemo.QQFieldElem)

Computes the orbit representatives of solution of the equation Qmat[v] = Xval.
If Xval = 0 then the orbit representatives of priitive isotropic vectors are computed.

# Examples

```
julia> eGram = Nemo.zero_matrix(Nemo.QQ, 3, 3); eGram[1,1] = 1; eGram[2,2] = -1; eGram[3,3] = -1;

julia> Indefinite.INDEF_FORM_GetOrbitRepresentative(eGram, QQ(0))
[1   0   -1]

```
"""
function INDEF_FORM_GetOrbitRepresentative(Qmat::Nemo.QQMatrix, Xval::Nemo.QQFieldElem)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  Xval_gap = GAP.evalstr(string("[", string(Xval), "]"))
  ListOrbitRepr_gap = GAP.Globals.INDEF_FORM_GetOrbitRepresentative(Qmat_gap, Xval_gap)
  return GAP.Globals.MatrixToOscar(ListOrbitRepr_gap)
end

"""
   INDEF_FORM_AutomorphismGroup(Qmat::Nemo.QQMatrix)

Computes a generating set of the automorphism group of the form Qmat

# Examples

```
julia> eGram = Nemo.matrix(Nemo.QQ, [0 1 0 0; 1 0 0 0; 0 0 -1 0; 0 0 0 -1])

julia> Indefinite.INDEF_FORM_AutomorphismGroup(eGram)
6-element Vector{QQMatrix}:
 [-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
 [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]
 [1 0 0 0; 1 1 1 1; -1 0 0 -1; -1 0 -1 0]
 [1 0 0 0; 2 1 2 0; -2 0 -1 0; 0 0 0 1]
 [1 1 1 1; 0 1 0 0; 0 -1 0 -1; 0 1 1 0]
 [5 4 6 2; 2 1 2 0; -4 -2 -4 -1; 2 2 3 0]

```
"""
function INDEF_FORM_AutomorphismGroup(Qmat::Nemo.QQMatrix)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  GRP_gap = GAP.Globals.INDEF_FORM_AutomorphismGroup(Qmat_gap)
  ListGens_gap = GAP.Globals.GeneratorsOfGroup(GRP_gap)
  return GAP.Globals.ListMatrixToOscar(ListGens_gap)
end

"""
   INDEF_FORM_TestEquivalence(Qmat1::Nemo.QQMatrix, Qmat2::Nemo.QQMatrix)

Test for equiavalence of two forms and if equivalent returns an equivalence

# Examples

```
julia> eGram1 = Nemo.matrix(Nemo.QQ, [0 1 0 0; 1 0 0 0; 0 0 -1 0; 0 0 0 -1])

julia> eGram2 = Nemo.matrix(Nemo.QQ, [2 1 0 0; 1 0 0 0; 0 0 -1 0; 0 0 0 -1])

julia> Indefinite.INDEF_FORM_TestEquivalence(eGram1, eGram2)
[1   1   0   0]
[1   0   0   0]
[0   0   1   0]
[0   0   0   1]

```
"""
function INDEF_FORM_TestEquivalence(Qmat1::Nemo.QQMatrix, Qmat2::Nemo.QQMatrix)
  Qmat1_gap = GAP.Globals.ReadOscarMatrix(Qmat1)
  Qmat2_gap = GAP.Globals.ReadOscarMatrix(Qmat2)
  TheEquiv_gap = GAP.Globals.INDEF_FORM_TestEquivalence(Qmat1_gap, Qmat2_gap)
  return GAP.Globals.MatrixToOscar(TheEquiv_gap)
end


"""
   INDEF_FORM_GetOrbit_IsotropicKplane(Qmat::Nemo.QQMatrix, k::Int64)

Computes the orbits of isotropic k-planes of the form Q

# Examples

```
julia> eGram = Nemo.matrix(Nemo.QQ, [0 1 0 0; 1 0 0 0; 0 0 0 2; 0 0 2 0])

julia> Indefinite.INDEF_FORM_GetOrbit_IsotropicKplane(eGram, 2)
????

```
"""
function INDEF_FORM_GetOrbit_IsotropicKplane(Qmat::Nemo.QQMatrix, k::Int64)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  ListRepr_gap = GAP.Globals.INDEF_FORM_GetOrbit_IsotropicKplane(Qmat_gap, k)
  return GAP.Globals.ListMatrixToOscar(ListRepr_gap)
end

"""
   INDEF_FORM_GetOrbit_IsotropicKplane(Qmat::Nemo.QQMatrix, k::Int64)

Computes the orbits of isotropic k-flags of the form Q

# Examples

```
julia> eGram = Nemo.matrix(Nemo.QQ, [0 1 0 0; 1 0 0 0; 0 0 0 2; 0 0 2 0])

julia> Indefinite.INDEF_FORM_GetOrbit_IsotropicKflag(eGram, 2)
????

```
"""
function INDEF_FORM_GetOrbit_IsotropicKflag(Qmat::Nemo.QQMatrix, k::Int64)
  Qmat_gap = GAP.Globals.ReadOscarMatrix(Qmat)
  ListRepr_gap = GAP.Globals.INDEF_FORM_GetOrbit_IsotropicKflag(Qmat_gap, k)
  return GAP.Globals.ListMatrixToOscar(ListRepr_gap)
end

