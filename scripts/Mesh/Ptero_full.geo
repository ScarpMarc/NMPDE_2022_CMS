// Gmsh project created on Tue Feb 14 23:24:55 2023
SetFactory("OpenCASCADE");

General.NumThreads = 0;
Geometry.OCCParallel = 1;
Geometry.OCCSewFaces = 1;

Merge "half_ptero_reducedx2.stp" ;//+

//HealShapes; // Just cause
/*DefineConstant[
    // Angle between two triangles above which an edge is considered as sharp
    angle = {40, Min 20, Max 120, Step 1,
      Name "Parameters/Angle for surface detection"},
    // For complex geometries, patches can be too complex, too elongated or too
    // large to be parametrized; setting the following option will force the
    // creation of patches that are amenable to reparametrization:
    forceParametrizablePatches = {0, Choices{0,1},
      Name "Parameters/Create surfaces guaranteed to be parametrizable"},
    // For open surfaces include the boundary edges in the classification process:
    includeBoundary = 1,
    // Force curves to be split on given angle:
    curveAngle = 180
  ];
  ClassifySurfaces{angle * Pi/180, includeBoundary, forceParametrizablePatches,
                   curveAngle * Pi / 180};
  
  // Create a geometry for all the discrete curves and surfaces in the mesh, by
  // computing a parametrization for each one
  CreateGeometry;*/

// Create ptero volume
Surface Loop(1) = Surface{:};
Volume(1) = {1};
//Physical Volume("PteroBody") = {1};
Physical Surface("PteroSurface") = {1:798}; // 798 for reduced mesh

BoundingBox;
bb() = BoundingBox Volume {1};

xmin = bb(0) * 1.1;//2;
ymin = bb(1) * 1.1; // 1.5;
zmin = bb(2) * 2;//12;
xmax = bb(3);
ymax = bb(4) * 1.2;// 7;
zmax = bb(5) * 1.5;//3.5;

Box(2) = {xmin,ymin,zmin, xmax-xmin,ymax-ymin,zmax-zmin};

Physical Surface("Inlet") = {801};
Physical Surface("Outlet") = {802};
Physical Surface("SlipSurface") = {800};
Physical Surface("Wall") = {799, 803, 804}; // 1 is the 
BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };
Physical Volume("SimulationField") = {3};

// Mesh

sampling_rate = 4;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.MinimumElementsPerTwoPi = 10;
Field[1] = Distance;
Field[1].SurfacesList = {1:798};
Field[1].Sampling = 2;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = sampling_rate / 80;
Field[2].SizeMax = sampling_rate;
Field[2].DistMin = .01;
Field[2].DistMax = .08;
Field[3] = MathEval;
Field[3].F = Sprintf("F1^3 + %g", sampling_rate / 150);
Field[7] = Min;
Field[7].FieldsList = {2, 3};
Background Field = 7;
/*
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
//Mesh.Algorithm = 5;
*/
Mesh.Algorithm3D = 10;

Mesh(3);

/*
Optimize the current mesh with the given algorithm (currently "Gmsh" for default tetrahedral mesh optimizer, "Netgen" for Netgen optimizer, "HighOrder" for direct high-order mesh optimizer, "HighOrderElastic" for high-order elastic smoother, "HighOrderFastCurving" for fast curving algorithm, "Laplace2D" for Laplace smoothing, "Relocate2D" and "Relocate3D" for node relocation).
*/
//OptimizeMesh "Gmsh";

Save "ptero.msh";

//Merge "half_ptero_reducedx2.stl" ;