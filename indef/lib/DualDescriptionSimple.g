__DualDescriptionLRS_Reduction:=function(EXT, GroupExt)
    local EXT_oscar, arith, method, method_oscar, ListIncd_oscar;
    EXT_oscar:=MatrixToOscar(EXT);
    method:="lrs_ring";
    method_oscar:=StringToOscar(method);
    ListIncd_oscar:=Julia.Indefinite.POLY_dual_description_group(method_oscar, EXT_oscar, GroupExt);
    return ReadOscarListIncd(ListIncd_oscar);
end;
