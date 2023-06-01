module Indefinite

import Oscar
import GAP
import Hecke
using Nemo
import polyhedral_jll
using BinaryWrappers

const indefinitedir = Base.pkgdir(Indefinite)
include("ExternalCalls.jl")
include("Functions.jl")


function __init__()
    binpaths = [ @generate_wrappers(polyhedral_jll) ]
    ENV["PATH"] = join([binpaths...,ENV["PATH"]], ":")
    GAP.Packages.load("grape")
    comm = string("\"", indefinitedir, "/indef/init.g\"")
    GAP.Globals.Read(GAP.evalstr(comm))
    comm = string("full_install_indefinite(\"", indefinitedir, "\")")
    GAP.evalstr(comm)
end

end
