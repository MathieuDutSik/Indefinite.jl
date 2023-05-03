using Nemo

#

FileSchlafli = "TestCases/GRP_LINPOLYTOPE_AUTOMORPHISM_case1_EXT27"
EXT27 = Indefinite.ReadMatrix_from_file(FileSchlafli)
GRP27 = Indefinite.GRP_LinPolytope_Automorphism(EXT27)

#

LATT_AUTO_case1_ListMat_E8 = Indefinite.ReadListMatrix_from_file("TestCases/LATT_AUTOMORPHISM_case1_ListMat_E8")
TheGRP1 = Indefinite.LATT_Automorphism(LATT_AUTO_case1_ListMat_E8)

#

LATT_ISOM_case1_ListM1 = Indefinite.ReadListMatrix_from_file("TestCases/LATT_ISOMORPHISM_case1_ListMat1_E8")
LATT_ISOM_case1_ListM2 = Indefinite.ReadListMatrix_from_file("TestCases/LATT_ISOMORPHISM_case1_ListMat2_2E8")
Equiv1 = Indefinite.LATT_Isomorphism(LATT_ISOM_case1_ListM1, LATT_ISOM_case1_ListM2)

LATT_ISOM_case2_ListM1 = Indefinite.ReadListMatrix_from_file("TestCases/LATT_ISOMORPHISM_case2_ListMat1_E8")
LATT_ISOM_case2_ListM2 = Indefinite.ReadListMatrix_from_file("TestCases/LATT_ISOMORPHISM_case2_ListMat2_E8_equiv")
Equiv2 = Indefinite.LATT_Isomorphism(LATT_ISOM_case2_ListM1, LATT_ISOM_case2_ListM2)

#

POLY_LP_case1_FAC  = Indefinite.ReadMatrix_from_file("TestCases/POLY_CDD_LINEARPROGRAMMING_case1_FAC")
POLY_LP_case1_Ineq = Indefinite.ReadVector_from_file("TestCases/POLY_CDD_LINEARPROGRAMMING_case1_Ineq")
ResultLP = Indefinite.LinearProgramming(POLY_LP_case1_FAC, POLY_LP_case1_Ineq)

#

POLY_DD_case1_EXT = Indefinite.ReadMatrix_from_file("TestCases/POLY_DUAL_DESCRIPTION_GROUP_case1_EXT")
POLY_DD_case1_GRP = Indefinite.ReadGroup_from_file("TestCases/POLY_DUAL_DESCRIPTION_GROUP_case1_GRP")
ListIncd = Indefinite.POLY_dual_description_group("lrs_ring", POLY_DD_case1_EXT, POLY_DD_case1_GRP)

#

POLY_SAMP_case1_EXT = Indefinite.ReadMatrix_from_file("TestCases/POLY_SAMPLING_FACETS_case1_FileI")
ListIncd = Indefinite.POLY_sampling_facets(POLY_SAMP_case1_EXT, "lp_cdd")

#

SV_NEAR_case1_Gram = Indefinite.ReadMatrix_from_file("TestCases/SV_NEAR_case1_Gram")
SV_NEAR_case1_V = Indefinite.ReadVector_from_file("TestCases/SV_NEAR_case1_V")
ListVect = Indefinite.LATT_near(SV_NEAR_case1_Gram, SV_NEAR_case1_V, Nemo.QQ(0))
