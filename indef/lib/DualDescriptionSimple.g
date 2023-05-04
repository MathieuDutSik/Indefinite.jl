__DualDescriptionLRS_Reduction:=function(EXT, GroupExt, ThePath)
    local EXT_oscar, GroupExt_oscar, arith, method, ListIncd_oscar;
    EXT_oscar:=MatrixToOscar(EXT);
    GroupExt_oscar:=PermutationGroupToOscar(Length(EXT), GroupExt);
    method:="lrs_ring";
    ListIncd_oscar:=Oscar.POLY_dual_description_group(method, EXT_oscar, GroupExt_oscar);
    return ReadOscarListIncd(ListIncd_oscar);
end;
