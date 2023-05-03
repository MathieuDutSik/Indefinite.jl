function parse_NN(strin::String)
  if length(strin) < 6
    val = parse(Int64, strin)
    return Nemo.QQ(val)
  else
    thepow = Nemo.QQ(1)
    len = length(strin)
    theret = Nemo.QQ(0)
    for i in 1:len
      echar = strin[len + 1 - i]
      val = Nemo.QQ(parse(Int64, echar))
      theret += val * thepow
      thepow = thepow * 10
    end
    return theret
  end
end


function parse_QQ(strin::String)
  if strin == "-"
    sign = -1
    strin_B = strin[2:end]
  else
    sign = 1
    strin_B = strin[1:end]
  end
  LStr = split(strin_B, "/")
  if length(LStr) == 1
    return sign * parse_NN(strin_B)
  else
    val1 = parse_NN(LStr[1])
    val2 = parse_NN(LStr[2])
    return sign * val1 / val2
  end
end

function WriteMatrix_to_stream(f::IOStream, M::Nemo.QQMatrix)
  n_rows = Nemo.nrows(M)
  n_cols = Nemo.ncols(M)
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

function WriteListMatrix_to_file(FileName::String, ListM::Vector{Nemo.QQMatrix})
  f = open(FileName, "w")
  WriteListMatrix_to_stream(f, ListM)
  close(f)
end

function ReadListMatrix_from_stream(f::IOStream)
  line = readline(f)
  n_mat = parse(Int64, line)
  ListM = Vector{Nemo.QQMatrix}(undef,0)
  for i_mat in 1:n_mat
    M = ReadMatrix_from_stream(f::IOStream)
    push!(ListM, M)
  end
  return ListM
end

function ReadListMatrix_from_file(FileName::String)
  f = open(FileName, "r")
  ListM = ReadListMatrix_from_stream(f::IOStream)
  close(f)
  return ListM
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
  line_first = readline(f)
  LStr = split(line_first, " ")
  n_row = parse(Int64, LStr[1])
  n_col = parse(Int64, LStr[2])
  print("n_row=", n_row, " n_col=", n_col, "\n")
  M = Nemo.zero_matrix(Nemo.QQ, n_row, n_col)
  for i in 1:n_row
    eline = readline(f)
    LStr = split(eline, " ")
    for j in 1:n_col
      val = parse_QQ(string(LStr[j+1]))
      M[i,j] = val
    end
  end
  return M
end

function ReadMatrix_from_file(FileName::String)
  f = open(FileName, "r")
  M = ReadMatrix_from_stream(f::IOStream)
  close(f)
  return M
end


function ReadVector_from_stream(f::IOStream)
  line_first = readline(f)
  n = parse(Int64, line_first)
  V = Nemo.zero_matrix(Nemo.QQ, 1, n)
  eline = readline(f)
  LStr = split(eline, " ")
  for j in 1:n
    val = parse_QQ(string(LStr[j+1]))
    M[1,j] = val
  end
  return V
end

function ReadVector_from_file(FileName::String)
  f = open(FileName, "r")
  V = ReadVector_from_stream(f::IOStream)
  close(f)
  return V
end

function ReadScalar_from_stream(f::IOStream)
  line_first = readline(f)
  return parse_QQ(line_first)
end



function WriteGroup_to_stream(f::IOStream, n, GRP::GAP.GAP_jll.GapObj)
  LGen = GAP.Globals.GeneratorsOfGroup(GRP)
  n_gen = GAP.Globals.Length(LGen)
  str_o = string(string(n), " ", string(n_gen), "\n")
  for i_gen in 1:n_gen
    eGen = LGen[i_gen]
    for i in 1:n
      eImg = GAP.Globals.OnPoints(i, eGen)
      str_o = string(str_o, " ", eImg - 1)
    end
    str_o = string(str_o, "\n")
  end
  write(f, str_o)
end

function WriteGroup_to_file(FileName::String, n, GRP::GAP.GAP_jll.GapObj)
  f = open(FileName, "w")
  WriteGroup_to_stream(f, n, GRP)
  close(f)
end


function ReadGroup_from_stream(f::IOStream)
  line_first = readline(f)
  LStr = split(line_first, " ")
  n = parse(Int64, LStr[1])
  n_gen = parse(Int64, LStr[2])
  print("n=", n, " n_gen=", n_gen, "\n")
  str_o = "Group("
  if n_gen > 0
    str_o = string(str_o, "[")
    for i_gen in 1:n_gen
      if i_gen > 1
        str_o = string(str_o, ",")
      end
      eline = readline(f)
      str_o = string(str_o, "PermList([")
      LStr = split(eline, " ")
      for i in 1:n
        if i > 1
          str_o = string(str_o, ",")
        end
        pos = parse(Int64, LStr[i+1])
        str_o = string(str_o, string(pos + 1))
      end
      str_o = string(str_o, "])")
    end
    str_o = string(str_o, "]")
  else
    str_o = string(str_o, "()")
  end
  str_o = string(str_o, ")")
  return GAP.evalstr(str_o)
end

function ReadGroup_from_file(FileName::String)
  f = open(FileName, "r")
  GRP = ReadGroup_from_stream(f::IOStream)
  close(f)
  return GRP
end


function GRP_LinPolytope_Automorphism(EXT::Nemo.QQMatrix)
  FileEXT = tempname()
  FileGroup = tempname()
  WriteMatrix_to_file(FileEXT, EXT)
  TheCommand = "GRP_LinPolytope_Automorphism"
  opt1 = "rational"
  opt2 = FileEXT
  opt3 = "Oscar"
  opt4 = FileGroup
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4`)
  GRP = ReadGroup_from_file(FileGroup)
  rm(FileEXT)
  rm(FileGroup)
  return GRP
end



function GRP_LinPolytope_Automorphism_GramMat(EXT::Nemo.QQMatrix, GramMat::Nemo.QQMatrix)
  FileEXT = tempname()
  FileGram = tempname()
  FileGroup = tempname()
  WriteMatrix_to_file(FileEXT, EXT)
  WriteMatrix_to_file(FileGram, GramMat)
  TheCommand = "GRP_LinPolytope_Automorphism_GramMat"
  opt1 = "rational"
  opt2 = FileEXT
  opt3 = FileGram
  opt4 = "Oscar"
  opt5 = FileGroup
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4 $opt5`)
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
  TheCommand = "GRP_LinPolytope_Isomorphism_GramMat"
  opt1 = FileEXT1
  opt2 = FileGram1
  opt3 = FileEXT2
  opt4 = FileGram2
  opt5 = "Oscar"
  opt6 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4 $opt5 $opt6`)
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
  TheCommand = "GRP_ListMat_Subset_EXT_Automorphism"
  opt1 = FileInput
  opt2 = "Oscar"
  opt3 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3`)
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
  TheCommand = "GRP_ListMat_Subset_EXT_Isomorphism"
  opt1 = FileInput
  opt2 = "Oscar"
  opt3 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3`)
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
  WriteMatrix_to_file(FileGram, GramMat)
  WriteVector_to_file(FileV, eV)
  TheCommand = "sv_near"
  opt1 = "rational"
  opt2 = choice
  opt3 = FileGram
  opt4 = FileV
  opt5 = "Oscar"
  opt6 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4 $opt5 $opt6`)
  MatVector = ReadMatrix_from_file(FileOut)
  rm(FileGram)
  rm(FileV)
  rm(FileOut)
  return MatVector
end

function POLY_dual_description_group(method::String, EXT::Nemo.QQMatrix, GRP::GAP.GAP_jll.GapObj)
  FileEXT = tempname()
  FileGRP = tempname()
  FileOut = tempname()
  WriteMatrix_to_file(FileEXT, EXT)
  n = nrows(EXT)
  WriteGroup_to_file(FileGRP, n, GRP)
  TheCommand = "POLY_dual_description_group"
  opt1 = "rational"
  opt2 = method
  opt3 = FileEXT
  opt4 = FileGRP
  opt5 = "Oscar"
  opt6 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4 $opt5 $opt6`)
  MatVector = ReadMatrix_from_file(FileOut)
  rm(FileEXT)
  rm(FileGRP)
  rm(FileOut)
  return MatVector
end

function IndefiniteReduction(GramMat::Nemo.QQMatrix)
  FileGram = tempname()
  FileOut = tempname()
  WriteMatrix_to_file(FileGram, GramMat)
  TheCommand = "IndefiniteReduction"
  opt1 = FileGram
  opt2 = "Oscar"
  opt3 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3`)
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
  TheCommand = "LATT_Automorphism"
  opt1 = FileListGram
  opt2 = "Oscar"
  opt3 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3`)
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
  TheCommand = "LATT_Isomorphism"
  opt1 = FileListGram1
  opt2 = FileListGram2
  opt3 = "Oscar"
  opt4 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4`)
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
  TheCommand = "POLY_cdd_LinearProgramming"
  opt1 = "rational"
  opt2 = FileFAC
  opt3 = FileIneq
  opt4 = "Oscar"
  opt5 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4 $opt5`)
  f = open(FileOut, "r")
  answer = readline(f)
  optimal_value = ReadScalar_from_stream(f)
  primal_solution = ReadVector_from_stream(f)
  dual_solution = ReadVector_from_stream(f)
  close(f)
  rm(FileFAC)
  rm(FileIneq)
  rm(FileOut)
  return [answer, optimal_value, primal_solution, dual_solution]
end

function POLY_samplingFacets(FAC::Nemo.QQMatrix, command::String)
  FileFAC = tempname()
  WriteMatrix_to_file(FileFAC, FAC)
  FileOut = tempname()
  TheCommand = "POLY_samplingFacets"
  opt1 = "rational"
  opt2 = command
  opt3 = FileFAC
  opt4 = "Oscar"
  opt5 = FileOut
  run(`$TheCommand $opt1 $opt2 $opt3 $opt4 $opt5`)
  ListIncd = ReadMatrix_from_file(FileOut)
  rm(FileFAC)
  rm(FileOut)
  return ListIncd
end

