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
  GRP = ReadGroupFile(FileGroup)
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

function GRP_ListMat_Subset_EXT_Automorphism(EXT::Nemo.QQMatrix, GramMat::Nemo.QQMatrix)
  FileEXT = tempname()
  FileGram = tempname()
  WriteMatrix_to_file(FileEXT, EXT)
  WriteMatrix_to_file(FileGram, GramMat)
  



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
  FileGram = tempname()
  FileOut = tempname()
  n = rows(GramMat)
  WriteMatrix_to_file(FileGram, GramMat)
  run(pipeline("IndefiniteReduction", FileGram, "Oscar", FileOut))
  B = zero_matrix(QQ, n, n)
  Mred = zero_matrix(QQ, n, n)
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