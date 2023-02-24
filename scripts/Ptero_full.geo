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

sampling_rate = 2;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.MinimumElementsPerTwoPi = 10;
Field[1] = Distance;
Field[1].SurfacesList = {1:798};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = sampling_rate / 50;
Field[2].SizeMax = sampling_rate;
Field[2].DistMin = 1;
Field[2].DistMax = 15;
Field[3] = MathEval;
Field[3].F = Sprintf("F1^3 + %g", sampling_rate / 60);
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
OptimizeMesh "Netgen";

Save "ptero.msh";

//Merge "half_ptero_reducedx2.stl" ;//+

/*
    For some strange reason Extrude does not properly create 
      a new volume: we will define it by hand by providing all 8 points
      instead of the first 4 only.
*/
/*
Point(100001) = {xmax, ymin, zmin, 1.0};
Point(100002) = {xmax, ymin, zmax, 1.0};
Point(100003) = {xmin, ymin, zmax, 1.0};
Point(100004) = {xmin, ymin, zmin, 1.0};
Line(100001) = {100001,100002}; Line(100002) = {100002,100003}; Line(100003) = {100003,100004}; Line(100004) = {100004,100001};
Curve Loop(100001) = {100001,100002,100003,100004};
Plane Surface(100002) = {100001}; // Surface 1 is the STL model
Physical Surface("Inlet") = {4};

// Extrude

Extrude {0, ymax-ymin, 0} {Surface{100002};}
Physical Surface("Outlet") = {5};
Physical Surface("SlipSurface") = {1};
Physical Surface("Wall") = {2, 3, 6}; // 1 is the body*/
/*

/*Point(5) = {xmax, ymax, zmin, 1.0};
Point(6) = {xmax, ymax, zmax, 1.0};
Point(7) = {xmin, ymax, zmax, 1.0};
Point(8) = {xmin, ymax, zmin, 1.0};
// Inlet
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(2) = {1}; // Surface 1 is the STL model
Physical Surface("Inlet") = {2};
// Outlet
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Curve Loop(2) = {5,6,7,8};
Plane Surface(3) = {2}; // Surface 1 is the STL model
Physical Surface("Outlet") = {3};
// Slip wall
Line(9) = {2,6}; Line(11) = {5,1};
// Other wall
Line(13) = {3,7}; Line(15) = {8,4};
// Create remaining 4 surfaces
Curve Loop(3) = {5, 9, 1, 11}; // Slip wall
Plane Surface(4) = {3}; // Surface 1 is the STL model
Physical Surface("SlipSurface") = {4};

Curve Loop(4) = {9,6,13,2}; // Upper surface
Plane Surface(5) = {4};
Curve Loop(5) = {7,15,3,13}; // Wall opposite to slip wall
Plane Surface(6) = {5}; 
Curve Loop(6) = {8,11,4,15}; // Lower surface
Plane Surface(7) = {6}; 
Physical Surface("Wall") = {1, 5, 6, 7}; // 1 is the body

Surface Loop(1) = {2, 3, 4, 5, 6, 7};
Volume(2) = {1};
BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };


//Surface Loop(2) = {2,3,4,5,6,7};
//Volume(2) = {2};


// Assign other names using automatic tags that
//  were created during extrusion


Physical Volume("SimulationField") = {3};
//Physical Volume("SimulationField") = {1};
*/
// Mesh
/*
sampling_rate = .1;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.MinimumElementsPerTwoPi = 10;
Field[1] = Distance;
Field[1].SurfacesList = {1};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = sampling_rate / 20;
Field[2].SizeMax = sampling_rate;
Field[2].DistMin = 1;
Field[2].DistMax = 15;
Field[3] = MathEval;
Field[3].F = Sprintf("F1^3 + %g", sampling_rate / 20);
Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 5;

// TODO: define volume for ptero body
Mesh(3);

/*
Optimize the current mesh with the given algorithm (currently "Gmsh" for default tetrahedral mesh optimizer, "Netgen" for Netgen optimizer, "HighOrder" for direct high-order mesh optimizer, "HighOrderElastic" for high-order elastic smoother, "HighOrderFastCurving" for fast curving algorithm, "Laplace2D" for Laplace smoothing, "Relocate2D" and "Relocate3D" for node relocation).
*/
// OptimizeMesh
//Save "ptero.msh";
