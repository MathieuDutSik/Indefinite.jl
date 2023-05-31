"""
   parse_NN(string::String)

Convert a string into a natural integer (negative values are excluded) implemented into Nemo.QQ

# Examples

```
julia> parse_NN("1234")
Nemo.QQ(1234)

```
"""
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


"""
   parse_QQ(string::String)

Convert a string into a rational number implemented into Nemo.QQ

# Examples

```
julia> parse_QQ("1234/23")
1234/23

julia> parse_QQ("-1234")
-1234

```
"""
function parse_QQ(strin::String)
  if strin[1] == '-'
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
    val1 = parse_NN(string(LStr[1]))
    val2 = parse_NN(string(LStr[2]))
    return sign * val1 / val2
  end
end

"""
   parse_QQ(string::String)

Convert a Nemo rational number into a string (without the "//" entry).
This is for the creation of input file for external programs.

# Examples

```
julia> parse_QQ("1234/23")
1234/23

julia> parse_QQ("-1234")
-1234

```
"""
function string_QQ(theval::Nemo.QQFieldElem)
  thestr = string(theval)
  LStr = split(thestr, "/")
  if length(LStr) == 1
    return thestr
  else
    return string(LStr[1], "/", LStr[3])
  end
end

"""
   WriteMatrix_to_stream(f::IOStream, M::Nemo.QQMatrix)

Write a matrix of rational numbers to a stream for use by external programs
"""
function WriteMatrix_to_stream(f::IOStream, M::Nemo.QQMatrix)
  n_rows = Nemo.nrows(M)
  n_cols = Nemo.ncols(M)
  str_o = string(string(n_rows), " ", string(n_cols), "\n")
  for i_row in 1:n_rows
    for i_col in 1:n_cols
      str_o = string(str_o, " ", string_QQ(M[i_row,i_col]))
    end
    str_o = string(str_o, "\n")
  end
  write(f, str_o)
end

"""
   WriteMatrix_to_file(FileName::String, M::Nemo.QQMatrix)

Write a matrix of rational numbers to a stream for use by external programs
"""
function WriteMatrix_to_file(FileName::String, M::Nemo.QQMatrix)
  f = open(FileName, "w")
  WriteMatrix_to_stream(f, M)
  close(f)
end

"""
   WriteListMatrix_to_stream(f::IOStream, ListM::Vector{Nemo.QQMatrix})

Write a vector of rational matrices to a stream for use by external programs
"""
function WriteListMatrix_to_stream(f::IOStream, ListM::Vector{Nemo.QQMatrix})
  n_mat = size(ListM)[1]
  write(f, string(n_mat), "\n")
  for i_mat in 1:n_mat
    WriteMatrix_to_stream(f, ListM[i_mat])
  end
end

"""
   WriteListMatrix_to_file(FileName::String, ListM::Vector{Nemo.QQMatrix})

Write a vector of rational matrices to a file for use by external programs
"""
function WriteListMatrix_to_file(FileName::String, ListM::Vector{Nemo.QQMatrix})
  f = open(FileName, "w")
  WriteListMatrix_to_stream(f, ListM)
  close(f)
end

"""
   ReadListMatrix_from_stream(f::IOStream)

Read a list of matrices from a stream
"""
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

"""
   ReadListMatrix_from_file(FileName::String)

Read a list of matrices from a stream
"""
function ReadListMatrix_from_file(FileName::String)
  f = open(FileName, "r")
  ListM = ReadListMatrix_from_stream(f::IOStream)
  close(f)
  return ListM
end

"""
   WriteVector_to_stream(f::IOStream, V::Nemo.QQMatrix)

Write a vector to a stream. The vector is coded as a matrix having
one row and N columns (with N the number of entries in the vector)
"""
function WriteVector_to_stream(f::IOStream, V::Nemo.QQMatrix)
  n_cols = Nemo.ncols(V)
  str_o = string(string(n_cols), "\n")
  for i_col in 1:n_cols
    str_o = string(str_o, " ", string_QQ(V[1,i_col]))
  end
  str_o = string(str_o, "\n")
  write(f, str_o)
end

"""
   WriteVector_to_file(FileName::String, V::Nemo.QQMatrix)

Write a vector to a file. The vector is coded as a matrix having
one row and N columns (with N the number of entries in the vector)
"""
function WriteVector_to_file(FileName::String, V::Nemo.QQMatrix)
  f = open(FileName, "w")
  WriteVector_to_stream(f::IOStream, V::Nemo.QQMatrix)
  close(f)
end

"""
   ReadMatrix_from_stream(f::IOStream)

Read a matrix from a stream. The entries should be in the following format:
nbRow nbCol
 A(1,1) .... A(1,nbCol)
.
.
 A(nbRow,1) .... A(nbRow,nbCol)
The spacing is important here. The lines of the matrix entry must start by a space
and each entry must be separated by only one space
"""
function ReadMatrix_from_stream(f::IOStream)
  line_first = readline(f)
  LStr = split(line_first, " ")
  n_row = parse(Int64, LStr[1])
  n_col = parse(Int64, LStr[2])
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

"""
   ReadMatrix_from_file(FileName::String)

Read a matrix from a file. The entries should be in the following format:
nbRow nbCol
 A(1,1) .... A(1,nbCol)
.
.
 A(nbRow,1) .... A(nbRow,nbCol)
The spacing is important here. The lines of the matrix entry must start by a space
and each entry must be separated by only one space
"""
function ReadMatrix_from_file(FileName::String)
  f = open(FileName, "r")
  M = ReadMatrix_from_stream(f::IOStream)
  close(f)
  return M
end

"""
   ReadVector_from_stream(f::IOStream)

Read a vector from a stream. The entries should be in the following format:
N
 A(1) .... A(N)
The spacing is important here. The lines of the vector entries must start by a space
and each entry must be separated by only one space
"""
function ReadVector_from_stream(f::IOStream)
  line_first = readline(f)
  n = parse(Int64, line_first)
  V = Nemo.zero_matrix(Nemo.QQ, 1, n)
  eline = readline(f)
  LStr = split(eline, " ")
  for j in 1:n
    val = parse_QQ(string(LStr[j+1]))
    V[1,j] = val
  end
  return V
end

"""
   ReadVector_from_file(FileName::String)

Read a vector from a file. The entries should be in the following format:
N
 A(1) .... A(N)
The spacing is important here. The lines of the vector entries must start by a space
and each entry must be separated by only one space
"""
function ReadVector_from_file(FileName::String)
  f = open(FileName, "r")
  V = ReadVector_from_stream(f::IOStream)
  close(f)
  return V
end

"""
   ReadScalar_from_stream(f::IOStream)

Read a single scalar from the stream. It mus occupy exactly one line of the stream.
"""
function ReadScalar_from_stream(f::IOStream)
  line_first = readline(f)
  return parse_QQ(line_first)
end



"""
   WriteGroup_to_stream(f::IOStream, n, GRP::GAP.GapObj)

Write a permutation group to a stream.
The integer n is the number of elements on which the group acts. This is needed
because GAP permutation groups do not have a specified number on elements on which
they act. The output is of the form
n n_gen
 0 1 2 ..... n-1
 1 0 2 ..... n-1
.
.
with n_gen the number of generators of the group. Each lines lists the images of
the elements 0, ...., n-1 by the group element. So the first line above corresponds
to the identity and the second line to the permutation (1,2) in GAP.
"""
function WriteGroup_to_stream(f::IOStream, n, GRP::GAP.GapObj)
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

"""
   WriteGroup_to_file(FileName::String, n, GRP::GAP.GapObj)

Write a permutation group to a file.
"""
function WriteGroup_to_file(FileName::String, n, GRP::GAP.GapObj)
  f = open(FileName, "w")
  WriteGroup_to_stream(f, n, GRP)
  close(f)
end

"""
   ReadGroup_from_stream(f::IOStream)

Write a permutation group from a stream.
"""
function ReadGroup_from_stream(f::IOStream)
  line_first = readline(f)
  LStr = split(line_first, " ")
  n = parse(Int64, LStr[1])
  n_gen = parse(Int64, LStr[2])
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

"""
   ReadGroup_from_stream(f::IOStream)

Write a permutation group from a file.
"""
function ReadGroup_from_file(FileName::String)
  f = open(FileName, "r")
  GRP = ReadGroup_from_stream(f::IOStream)
  close(f)
  return GRP
end

"""
   GRP_LinPolytope_Automorphism(EXT::Nemo.QQMatrix)

Computes the group of linear automorphism preserving the matrix put as input.
"""
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

"""
   GRP_LinPolytope_Automorphism_GramMat(EXT::Nemo.QQMatrix, GramMat::Nemo.QQMatrix)

Computes the group of linear automorphism preserving the matrix EXT and the Gram matrix put on input.
"""
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

"""
   GRP_LinPolytope_Isomorphism_GramMat(EXT1::Nemo.QQMatrix, GramMat1::Nemo.QQMatrix, EXT2::Nemo.QQMatrix, GramMat2::Nemo.QQMatrix)

Computes the equivalence of EXT1 with EXT2 mapping GramMat1 to GramMat2 if existing. The returned vector maps the vertices.
If its length is zero then there is no equivalence
"""
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
  TheEquiv = ReadVector_from_file(FileOut)
  rm(FileEXT1)
  rm(FileGram1)
  rm(FileEXT2)
  rm(FileGram2)
  rm(FileOut)
  return TheEquiv
end

"""
   GRP_ListMat_Subset_EXT_Automorphism(EXT::Nemo.QQMatrix, ListGramMat::Vector{Nemo.QQMatrix}, Vdiag::Nemo.QQMatrix)

Computes the linear automorphism group of EXT preserving the list of Gram matrices and the vertex weight in Vdiag.
"""
function GRP_ListMat_Subset_EXT_Automorphism(EXT::Nemo.QQMatrix, ListGramMat::Vector{Nemo.QQMatrix}, Vdiag::Nemo.QQMatrix)
  FileInput = tempname()
  FileGroup = tempname()
  f = open(FileInput, "w")
  WriteListMatrix_to_stream(f, ListGramMat)
  WriteMatrix_to_stream(f, EXT)
  WriteVector_to_stream(f, Vdiag)
  close(f)
  TheCommand = "GRP_ListMat_Subset_EXT_Automorphism"
  opt1 = FileInput
  opt2 = "Oscar"
  opt3 = FileGroup
  run(`$TheCommand $opt1 $opt2 $opt3`)
  GRP = ReadGroup_from_file(FileGroup)
  rm(FileInput)
  rm(FileGroup)
  return GRP
end

"""
   GRP_ListMat_Subset_EXT_Invariant(EXT::Nemo.QQMatrix, ListGramMat::Vector{Nemo.QQMatrix}, Vdiag::Nemo.QQMatrix)

Computes an invariant under the equivalence preserving list of gram matrices and vertex weight
"""
function GRP_ListMat_Subset_EXT_Invariant(EXT::Nemo.QQMatrix, ListGramMat::Vector{Nemo.QQMatrix}, Vdiag::Nemo.QQMatrix)
  FileInput = tempname()
  FileOut = tempname()
  f = open(FileInput, "w")
  WriteListMatrix_to_stream(f, ListGramMat)
  WriteMatrix_to_stream(f, EXT)
  WriteVector_to_stream(f, Vdiag)
  close(f)
  TheCommand = "GRP_ListMat_Subset_EXT_Invariant"
  opt1 = FileInput
  opt2 = FileOut
  run(`$TheCommand $opt1 $opt2`)
  f = open(FileOut, "r")
  result = readline(f)
  close(f)
  rm(FileInput)
  rm(FileOut)
  return result
end

"""
   GRP_ListMat_Subset_EXT_Isomorphism(EXT1::Nemo.QQMatrix, ListGramMat1::Vector{Nemo.QQMatrix}, Vdiag1::Nemo.QQMatrix, EXT2::Nemo.QQMatrix, ListGramMat2::Vector{Nemo.QQMatrix}, Vdiag2::Nemo.QQMatrix)

Test for linear isomorphism of EXT1 and EXT2 mapping ListGramMat1 to ListGramMat2 and mapping the vertex weight Vdiag1 to Vdiag2
"""
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

"""
   LATT_near(GramMat::Nemo.QQMatrix, eV::Nemo.QQMatrix, Dist_mat::Nemo.QQMatrix)

Computes the vertors in the lattice up to distance Dist. If Dist=0 then the
shortest vectors are returned.
"""
function LATT_near(GramMat::Nemo.QQMatrix, eV::Nemo.QQMatrix, Dist_mat::Nemo.QQMatrix)
  Dist = Dist_mat[1,1]
  if Dist == 0
    choice = "nearest"
  else
    choice = string("near=", string_QQ(Dist))
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

"""
   POLY_dual_description_group(method::String, EXT::Nemo.QQMatrix, GRP::GAP.GapObj)

Computes the orbits of EXT for the group GRP. The method used is given in method
and can be lrs_ring or cdd.
"""
function POLY_dual_description_group(method::String, EXT::Nemo.QQMatrix, GRP::GAP.GapObj)
  FileEXT = tempname()
  FileGRP = tempname()
  FileOut = tempname()
  WriteMatrix_to_file(FileEXT, EXT)
  n = Nemo.nrows(EXT)
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

"""
   IndefiniteReduction(GramMat::Nemo.QQMatrix)

Compute a reduce form of the indefinite matrix GramMat.
It does not have to be positive definite and the result
is not canonical in any way. But it should be fast and
the result should have much smaller coefficient than the
original matrix.
"""
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

"""
   LATT_Automorphism(ListGramMat::Vector{Nemo.QQMatrix})

Compute the list of matrix generators of the group of matrices
preserving the list of Gram matrices. The first matrix has to
be positive definite.
"""
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

"""
   LATT_Isomorphism(ListGramMat1::Vector{Nemo.QQMatrix}, ListGramMat2::Vector{Nemo.QQMatrix})

Computes the integral equivalence between ListGramMat1 and ListGramMat2 if one exists. If one
exist then an integral matrix is returned, if not then the returned matrix has 0 lines and rows.
"""
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

"""
   LinearProgramming(InequalitySet::Nemo.QQMatrix, ToBeMinimized::Nemo.QQMatrix)

Computes the Linear programming for the set of inequalities and the minimization
of ToBeMinimized.
"""
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

"""
   POLY_sampling_facets(FAC::Nemo.QQMatrix, command::String)

Sampling vertices of the polytope defined by the inequaltiies of FAC.
The command used is described in the input variable command.
"""
function POLY_sampling_facets(FAC::Nemo.QQMatrix, command::String)
  FileFAC = tempname()
  WriteMatrix_to_file(FileFAC, FAC)
  FileOut = tempname()
  TheCommand = "POLY_sampling_facets"
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

