# Run Elder-McCammond's cat.g on the two shapes used in the census.
Read("cat.g");
verbose := false;
Print("=== disphenoid (coxeter shape [4,3,3,3,3,4]) ===\n");
coxeter();
test();
Print("DISPHENOID vGalleries: ", Length(vGalleries), "  mGalleries: ", Length(mGalleries), "\n");
Print("vGalleries := ", vGalleries, ";\n");
Print("mGalleries := ", mGalleries, ";\n");
Print("=== cube orthoscheme [1,2,3,1,2,1] ===\n");
createTetra([[1,2,3,1,2,1]]);
test();
Print("CUBE vGalleries: ", Length(vGalleries), "  mGalleries: ", Length(mGalleries), "\n");
Print("vGalleries := ", vGalleries, ";\n");
Print("mGalleries := ", mGalleries, ";\n");
QUIT;
