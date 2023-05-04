module Indefinite

import Oscar
import GAP
import Hecke
import Nemo
import polyhedral_jll
using BinaryWrappers

include("ExternalCalls.jl")


function __init__()
    binpaths = [ @generate_wrappers(polyhedral_jll) ]
    ENV["PATH"] = join([binpaths...,ENV["PATH"]], ":")
    GAP.Packages.load("grape")
    GAP.Globals.Read(GAP.evalstr("\"indef/init.g\""))
end

end
