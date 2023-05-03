module Indefinite

import GAP
import Hecke
import Nemo
import polyhedral_jll
using BinaryWrappers

include("ExternalCalls.jl")


function __init__()
    binpaths = [ @generate_wrappers(polyhedral_jll) ]
    ENV["PATH"] = join([binpaths...,ENV["PATH"]], ":")
end

end
