FileSchlafli = "TestCases/GRP_LINPOLYTOPE_AUTOMORPHISM_case1_EXT27"

EXT27 = Indefinite.ReadMatrix_from_file(FileSchlafli)

GRP27 = Indefinite.GRP_LinPolytope_Automorphism(EXT27)

#

LATT_AUTO_case1_ListMat_E8 = Indefinite.ReadListMatrix_from_file("TestCases/LATT_AUTOMORPHISM_case1_ListMat_E8")
TheGRP1 = Indefinite.LATT_Automorphism(LATT_AUTO_case1_ListMat_E8)