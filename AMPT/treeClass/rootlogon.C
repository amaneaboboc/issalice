{
gInterpreter->AddIncludePath("/alice/treeClass");
gSystem->AddIncludePath("-I. -I /alice/treeClass");
gSystem->AddDynamicPath(".:/alice/treeClass");

gSystem->Load("/alice/treeClass/libMyTree");
//gSystem->ListLibraries();
}
