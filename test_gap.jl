using Nemo
using GAP

#

CaseRepr1 = false
if CaseRepr1
   eGram = Indefinite.ReadMatrix_from_file("TestLor/U_I3_mat1")
   ListRepr = Indefinite.INDEF_FORM_GetOrbitRepresentative(eGram, QQ(0))
end

#

CaseAuto1 = false
if CaseAuto1
   eGram = Indefinite.ReadMatrix_from_file("TestLor/Gram_U+I3")
   ListGens = Indefinite.INDEF_FORM_AutomorphismGroup(eGram)
end

CaseAuto2 = false
if CaseAuto2
   eGram = Indefinite.ReadMatrix_from_file("TestLor/U_E8_mat1")
   ListGens = Indefinite.INDEF_FORM_AutomorphismGroup(eGram)
end

CaseAuto3 = false
if CaseAuto3
   eGram = Indefinite.ReadMatrix_from_file("TestLor/U_2U_2I3_mat1")
   ListGens = Indefinite.INDEF_FORM_AutomorphismGroup(eGram)
end

#

CaseEquiv1 = false
if CaseEquiv1
   Gram1 = Indefinite.ReadMatrix_from_file("TestLor/U_I3_mat1")
   Gram2 = Indefinite.ReadMatrix_from_file("TestLor/U_I3_mat2")
   TheEquiv = Indefinite.INDEF_FORM_TestEquivalence(Gram1, Gram2)
end

CaseEquiv2 = false
if CaseEquiv2
   Gram1 = Indefinite.ReadMatrix_from_file("TestLor/U_E8_mat1")
   Gram2 = Indefinite.ReadMatrix_from_file("TestLor/U_E8_mat2")
   TheEquiv = Indefinite.INDEF_FORM_TestEquivalence(Gram1, Gram2)
end

CaseEquiv3 = false
if CaseEquiv3
   Gram1 = Indefinite.ReadMatrix_from_file("TestLor/U_2U_2I3_mat1")
   Gram2 = Indefinite.ReadMatrix_from_file("TestLor/U_2U_2I3_mat2")
   TheEquiv = Indefinite.INDEF_FORM_TestEquivalence(Gram1, Gram2)
end

#

CaseIsoPlane1 = false
if CaseIsoPlane1
   Gram1 = Indefinite.ReadMatrix_from_file("TestLor/U_2U_2I3_mat1")
   ListRepr = Indefinite.INDEF_FORM_GetOrbit_IsotropicKplane(Gram1, 2)
end

#

CaseIsoFlag1 = true
if CaseIsoFlag1
   Gram1 = Indefinite.ReadMatrix_from_file("TestLor/U_2U_2I3_mat1")
   ListRepr = Indefinite.INDEF_FORM_GetOrbit_IsotropicKflag(Gram1, 2)
end
