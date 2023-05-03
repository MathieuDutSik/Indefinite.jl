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

