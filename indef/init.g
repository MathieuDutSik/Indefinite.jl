# DeclareAutoPackage("indefinite", "1.0", true);
# DeclarePackageAutoDocumentation("indefinite", "doc");


full_install_indefinite:=function(indefinitedir)
    local FCT;
    FCT:=function(filename)
        return Concatenation(indefinitedir, "/indef/lib/", filename);
    end;


    Read(FCT("InputOutput.g"));
    Read(FCT("Fractions.g"));
    Read(FCT("Fundamental.g"));
    Read(FCT("PositivitySymmetricMatrices.g"));
    Read(FCT("GroupAction.g"));
    Read(FCT("MyGraphicalFunctions.g"));
    Read(FCT("AutomorphismPolytope.g"));
    Read(FCT("LinearProgramming.g"));
    Read(FCT("SetFunctionality.g"));
    Read(FCT("DualDescriptionSimple.g"));
    Read(FCT("DualDescriptionAndADM.g"));
    Read(FCT("CVP.g"));
    Read(FCT("LatticeIsomorphy.g"));
    Read(FCT("LatticeDelaunays.g"));
    Read(FCT("IndefiniteFormsFundamental.g"));
    Read(FCT("Lorentzian.g"));
    Read(FCT("IndefiniteForms.g"));
end;
