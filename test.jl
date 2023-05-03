FileSchlafli = "TestCases/GRP_LINPOLYTOPE_AUTOMORPHISM_case1_EXT27"

EXT27 = Indefinite.ReadMatrix_from_file(FileSchlafli)

GRP27 = Indefinite.GRP_LinPolytope_Automorphism(EXT27)
