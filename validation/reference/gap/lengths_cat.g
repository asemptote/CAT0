# Geodesic length data for each vertex gallery found by cat.g (cube orthoscheme).
Read("cat.g");
verbose := false;
createTetra([[1,2,3,1,2,1]]);
test();
Print("tList := ", tList, ";\n");
for gg in vGalleries do
  deeplyConsider(gg);
  j := dotProd(gCoords[1], gCoords[gLength+2]);
  a := tList[AbsInt(g[1])][1];  b := tList[g[gLength]][3];
  Print("gallery ", gg, "  tetrahedra=", Length(gg), "  start/end sq-lengths ", [a, b],
        "  dot=", algToText(j), "  cos^2=", algToNum(algSMult(1/(a*b), algMult(j,j))),
        "  sign(dot)=", signNum(j), "\n");
od;
QUIT;
