

FuncRandomDirection:=function(n, siz)
  local eList, i;
  eList:=[0];
  for i in [1..n]
  do
    Add(eList, Random([-siz..siz]));
  od;
  return eList;
end;
