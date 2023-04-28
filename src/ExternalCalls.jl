function WriteMatrix_to_stream(f::IOStream, M::Nemo.QQMatrix)
  n_rows = rows(M)
  n_cols = cols(M)
  str_o = string(string(n_rows), " ", string(n_cols), "\n")
  for i_row in 1:n_rows
    for i_col in 1:n_cols
      str_o = string(str_o, " ", string(M[i_row,i_col]))
    end
    str_o = string(str_o, "\n")
  end
  write(f, str_o)
end

function WriteMatrix_to_file(FileName::String, M::Nemo.QQMatrix)
  f = open(FileName, "w")
  WriteMatrix_to_stream(f, M)
  close(f)
end

function WriteListMatrix_to_stream(f::IOStream, ListM::Vector{Nemo.QQMatrix})
  n_mat = size(ListM)[1]
  write(f, string(n_mat))
  for i_mat in 1:n_mat
    WriteMatrix_to_stream(f, ListM[i_mat])
  end
end

function WriteVector_to_stream(f::IOStream, V::Nemo.QQMatrix)
  n_cols = cols(M)
  str_o = string(string(n_cols), "\n")
  for i_col in 1:n_cols
    str_o = string(str_o, " ", string(M[1,i_col]))
  end
  str_o = string(str_o, "\n")
  write(f, str_o)
end

function WriteVector_to_file(FileName::String, V::Nemo.QQMatrix)
  f = open(FileName, "w")
  WriteVector_to_stream(f::IOStream, V::Nemo.QQMatrix)
  close(f)
end





function ReadMatrix_from_stream(f::IOStream)
  n_row = read(f, Int64)
  n_col = read(f, Int64)
  M = zero_matrix(QQ, n_row, n_col)
  for i in 1:n_row
    for j in 1:n_col
      M[i,j] = read(f, QQ)
    end
  end
  return M
end


function ReadGroup_from_stream(f::IOStream)
  n = read(f, Int64)
  n_gen = read(f, Int64)
  str_o = "Group("
  if n_gen > 0
    str_o = string(str_o, "[")
    for i_gen in 1::n_gen
    do
      if i_gen > 1
        str_o = string(str_o, ",")
      end
      str_o = string(str_o, "PermList([")
      for i in 1:n
        pos = read(f, Int64)
        if i > 1
          str_o = string(str_o, ",")
        end
        str_o = string(str_o, string(pos + 1)))
      end
      str_o = string(str_o, "])")
    od
    str_o = string(str_o, "]")
  else
    str_o = string(str_o, "()")
  end
  str_o = string(str_o, ")")
  return GAP.evalstr(str_o)
end

function ReadGroup_from_file(FileName::String)
  f = open(FileOut, "r")
  GRP = ReadGroup_from_stream(f::IOStream)
  close(f)
  return GRP
end


function GRP_LinPolytope_Automorphism_GramMat(EXT::Nemo.QQMatrix, GramMat::Nemo.QQMatrix)
  FileEXT = tempname()
  FileGram = tempname()
  FileGroup = tempname()
  WriteMatrix_to_file(FileEXT, EXT)
  WriteMatrix_to_file(FileGram, GramMat)
  run(pipeline("GRP_LinPolytope_Automorphism_GramMat", "rational", FileEXT, FileGram, FileGroup))
  GRP = ReadGroup_from_file(FileGroup)
  rm(FileEXT)
  rm(FileGram)
  rm(FileGroup)
  return GRP
end

function GRP_LinPolytope_Isomorphism_GramMat(EXT1::Nemo.QQMatrix, GramMat1::Nemo.QQMatrix, EXT2::Nemo.QQMatrix, GramMat2::Nemo.QQMatrix)
  FileEXT1 = tempname()
  FileGram1 = tempname()
  FileEXT2 = tempname()
  FileGram2 = tempname()
  FileOut = tempname()
  WriteMatrix_to_file(FileEXT1, EXT1)
  WriteMatrix_to_file(FileGram1, GramMat1)
  WriteMatrix_to_file(FileEXT2, EXT2)
  WriteMatrix_to_file(FileGram2, GramMat2)
  run(pipeline("GRP_LinPolytope_Isomorphism_GramMat", FileEXT1, FileGram1, FileEXT2, FileGram2, "Oscar", FileOut))
  TheEquiv = ReadMatrix_from_file(FileOut)
  rm(FileEXT1)
  rm(FileGram1)
  rm(FileEXT2)
  rm(FileGram2)
  rm(FileOut)
  return TheEquiv
end

function GRP_ListMat_Subset_EXT_Automorphism(EXT::Nemo.QQMatrix, ListGramMat::Vector{Nemo.QQMatrix}, Vdiag::Nemo.QQMatrix)
  FileInput = tempname()
  FileOut = tempname()
  f = open(FileInput, "w")
  WriteListMatrix_to_stream(f, ListGramMat)
  WriteMatrix_to_stream(f, EXT)
  WriteVector_to_stream(f, Vdiag)
  close(f)
  run(pipeline("GRP_ListMat_Subset_EXT_Automorphism", FileInput, "Oscar", FileOut))
  GRP = ReadGroup_from_file(FileGroup)
  rm(FileInput)
  rm(FileOut)
  return GRP
end

function GRP_ListMat_Subset_EXT_Isomorphism(EXT1::Nemo.QQMatrix, ListGramMat1::Vector{Nemo.QQMatrix}, Vdiag1::Nemo.QQMatrix, EXT2::Nemo.QQMatrix, ListGramMat2::Vector{Nemo.QQMatrix}, Vdiag2::Nemo.QQMatrix)
  FileInput = tempname()
  FileOut = tempname()
  f = open(FileInput, "w")
  WriteListMatrix_to_stream(f, ListGramMat1)
  WriteMatrix_to_stream(f, EXT1)
  WriteVector_to_stream(f, Vdiag1)
  WriteListMatrix_to_stream(f, ListGramMat2)
  WriteMatrix_to_stream(f, EXT2)
  WriteVector_to_stream(f, Vdiag2)
  close(f)
  run(pipeline("GRP_ListMat_Subset_EXT_Isomorphism", FileInput, "Oscar", FileOut))
  eVect = ReadVector_from_file(FileOut)
  rm(FileInput)
  rm(FileOut)
  return eVect
end

function LATT_near(GramMat::Nemo.QQMatrix, eV::Nemo.QQMatrix, Dist::Nemo.QQFieldElem)
  if Dist == 0
    choice = "nearest"
  else
    choice = string("near=", string(Dist))
  end
  FileGram = tempname()
  FileV = tempname()
  FileOut = tempname()
  WriteMatrix_from_file(FileGram, GramMat)
  WriteMatrix_from_file(FileV, eV)
  run(pipeline("LATT_near", "rational", choice, FileGram, FileV, "Oscar", FileOut))
  MatVector = ReadMatrix_from_file(FileOut)
  rm(FileGram)
  rm(FileV)
  rm(FileOut)
  return MatVector
end

function POLY_dual_description_group(method::String, EXT::Nemo.QQMatrix, GroupExt::Nemo.ZZMatrix)
  return "not done"
end

function IndefiniteReduction(GramMat::Nemo.QQMatrix)
  FileGram = tempname()
  FileOut = tempname()
  WriteMatrix_to_file(FileGram, GramMat)
  run(pipeline("IndefiniteReduction", FileGram, "Oscar", FileOut))
  #
  f = open(FileOut, "r")
  B = ReadMatrix_from_stream(f)
  Mred = ReadMatrix_from_stream(f)
  close(f)
  rm(FileGram)
  rm(FileOut)
  return [B, Mred]
end

function LATT_Automorphism(ListGramMat::Vector{Nemo.QQMatrix})
  FileListGram = tempname()
  FileOut = tempname()
  WriteListMatrix_to_file(FileListGram, ListGramMat)
  run(pipeline("LATT_Automorphism", FileListGram, "Oscar", FileOut))
  ListGens = ReadListMatrix_from_file(FileOut)
  rm(FileListGram)
  rm(FileOut)
  return ListGens
end

function LATT_Isomorphism(ListGramMat1::Vector{Nemo.QQMatrix}, ListGramMat2::Vector{Nemo.QQMatrix})
  FileListGram1 = tempname()
  FileListGram2 = tempname()
  FileOut = tempname()
  WriteListMatrix_to_file(FileListGram1, ListGramMat1)
  WriteListMatrix_to_file(FileListGram2, ListGramMat2)
  run(pipeline("LATT_Isomorphism", FileListGram1, FileListGram2, "Oscar", FileOut))
  TheEquiv = ReadMatrix_from_file(FileOut)
  rm(FileListGram1)
  rm(FileListGram2)
  rm(FileOut)
  return TheEquiv
end

function LinearProgramming(InequalitySet::Nemo.QQMatrix, ToBeMinimized::Nemo.QQMatrix)
  FileFAC = tempname()
  FileIneq = tempname()
  WriteMatrix_to_file(FileFAC, InequalitySet)
  WriteVector_to_file(FileIneq, ToBeMinimized)
  FileOut = tempname()
  run(pipeline("POLY_cdd_LinearProgramming", "rational", FileFAC, FileIneq, "Oscar", FileOut))
  f = open(FileOut, "r")
  answer = read(f, String)
  optimal_value = read(f, QQFieldElem)
  primal_solution = ReadVector_from_stream(f)
  dual_solution = ReadVector_from_stream(f)
  close(f)
  rm(FileFAC)
  rm(FileIneq)
  rm(FileOut)
  return [optimal_value, primal_solution, dual_solution]
end

function POLY_samplingFacets(FAC::Nemo.QQMatrix, command::String)
  FileFAC = tempname()
  WriteMatrix_to_file(FileFAC, FAC)
  FileOut = tempname()
  run(pipeline("POLY_samplingFacets", "rational", command, FileFAC, "Oscar", FileOut))
  ListIncd = ReadMatrix_from_file(FileOut)
  rm(FileFAC)
  rm(FileOut)
  return ListIncd
end