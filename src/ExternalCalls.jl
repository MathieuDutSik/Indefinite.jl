function WriteMatrix_to_file(FileName::String, M::Nemo.QQMatrix)
  n_rows = rows(M)
  n_cols = cols(M)
  str_o = string(string(n_rows), " ", string(n_cols), "\n")
  for i_row in 1:n_rows
    for i_col in 1:n_cols
      str_o = string(str_o, " ", string(M[i_row,i_col]))
    end
    str_o = string(str_o, "\n")
  end
  f = open(FileName, "w")
  write(f, str_o)
  close(f)
end



function GRP_LinPolytope_Automorphism_GramMat(EXT::Nemo.QQMatrix, GramMat::Nemo.QQMatrix)
  FileEXT = tempname()
  FileGram = tempname()
  WriteMatrix_to_file(FileEXT, EXT)
  WriteMatrix_to_file(FileGram, GramMat)
  
  return "not done"
end

function GRP_LinPolytope_Isomorphism_GramMat(EXT1::Nemo.QQMatrix, GramMat1::Nemo.QQMatrix, EXT2::Nemo.QQMatrix, GramMat2::Nemo.QQMatrix)
  return "not done"
end

function GRP_ListMat_Subset_EXT_Automorphism(EXT::Nemo.QQMatrix, GramMat::Nemo.QQMatrix)
  return "not done"
end

function GRP_ListMat_Subset_EXT_Isomorphism(EXT1::Nemo.QQMatrix, GramMat1::Nemo.QQMatrix, EXT2::Nemo.QQMatrix, GramMat2::Nemo.QQMatrix)
  return "not done"
end

function LATT_near(GramMat::Nemo.QQMatrix, eV::Nemo.QQMatrix, Dist::Nemo.QQFieldElem)
  return "not done"
end

function POLY_dual_description_group(method::String, EXT::Nemo.QQMatrix, GroupExt::Nemo.ZZMatrix)
  return "not done"
end

function IndefiniteReduction(GramMat::Nemo.QQMatrix)
  return "not done"
end

function LATT_Automorphism(ListGramMat::Vector{Nemo.QQMatrix})
  return "notdone"
end

function LATT_Isomorphism(ListGramMat1::Vector{Nemo.QQMatrix}, ListGramMat2::Vector{Nemo.QQMatrix})
  return "not done"
end

function LinearProgramming(InequalitySet::Nemo.QQMatrix, ToBeMinimized::Nemo.QQMatrix)
  return "not done"
end

function POLY_samplingFacets(FAC::Nemo.QQMatrix, command::String)
  return "not done"
end