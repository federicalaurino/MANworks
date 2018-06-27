//Geometrical parameters

// Define the coordinate range of the mesh
x0 = -0.25;
x1 = +0.25;
y0 = -1.0;
y1 = +1.0;
z0 = -1.0;
z1 = +1.0;

// Define the center of the vessel (directed along x-axis)
yc= 0.0;
zc= 0.0;

// Radius of the vessel
R = 0.05;
// number of elements for unit size
Num_el= 32.0;
// mesh element size at each point
point_size = 0.1;


//Define the 3d tissue slab

Nblayers=  Num_el * (y1-y0) ;
Point(1) = {x0, y0, z0, point_size};
Extrude {0, y1-y0, 0} {
  Point{1}; Layers{ Nblayers };
}
//Point(2) = {x0, y1, z0, point_size};
//Line(1) = {1, 2};

Nblayers=  Num_el * (z1-z0) ;
Extrude {0, 0, z1-z0} {
  Line{1}; Layers{ Nblayers };
}

Nblayers=Num_el* (x1-x0) ;
Extrude {x1-x0, 0, 0} {
  Surface{5}; Layers {Nblayers};
}

// define the vessel
Point(101) = {x0, yc, zc, point_size};
Point(102) = {x1, yc, zc, point_size};

Point(103) = {x0, yc+R, zc, point_size};
Point(104) = {x1, yc+R, zc, point_size};

Nblayers=Num_el* R ;
Extrude {0, R, 0} {
  Point{101}; Layers{ Nblayers };
}


Nblayers=Num_el* Pi*R/2 ;
Extrude {{1, 0, 0}, {x0, yc, zc}, Pi/2} {
  Line{28};Layers {Nblayers};
}
Extrude {{1, 0, 0}, {x0, yc, zc}, Pi/2} {
  Line{29};Layers {Nblayers};
}
Extrude {{1, 0, 0}, {x0, yc, zc}, Pi/2} {
  Line{32};Layers {Nblayers};
}
Extrude {{1, 0, 0}, {x0, yc, zc}, +Pi/2} {
  Line{35};Layers {Nblayers};
}


Nblayers=Num_el* (x1-x0) ;
Extrude {x1-x0, 0, 0} {
  Surface{31, 34, 37, 40};Layers {Nblayers};
}

// Define Physical groups for boundaries 

// NEUMANN CONDITION FACES (trasversal faces to the vessel)
//Physical Surface(184) = {27, 5};
// DIRICHLET CONDITION FACES (longitudinal faces to the vessel)
//Physical Surface(185) = {22, 18, 14, 26};
// Whole VOLUME
Physical Volume(186) = {1, 3, 4, 5, 2};
