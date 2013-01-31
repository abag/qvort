// To make the figure use the following (using desired width (W) and height (H))
// povray +IplotdecompoBox.pov +OtestBox.png +FN +W800 +H600 +V -D +X
// To remove the background add option +UA


////global_settings { assumed_gamma 2.2 } DO NOT USE
global_settings { assumed_gamma 1.0 }
global_settings { max_trace_level 30 }
global_settings { ambient_light rgb <1, 1, 1> }

#include "colors.inc"    // The include files contain
#include "stones.inc"    // pre-defined scene elements
#include "textures.inc"  
#include "shapes.inc"
#include "glass.inc"
#include "metals.inc"
#include "woods.inc"

#declare iorCylinder = 1.5;   // ior value for the cylinder glass 1.5
#declare DD = 0.5;  // Thickness for the walls
#declare LL = 10.08; //Box Size
#declare LL2 = LL+DD;

// First change the left-handed coordinate system to right-handed by commands up, right, and sky:
camera {
location <4*LL,3*LL,2*LL>
up <0,1,0>           // To get right-handed coordinate system
right <-1.33,0,0>    // To get right-handed coordinate system
// right <-1.00,0,0> // To get right-handed coordinate system
look_at <0,0,0>
sky <0,0,1>          // To get right-handed coordinate system
angle 24
}


// background { color White }
background { color Black }

light_source { <3*LL, 0*LL, 2*LL> color red 0.6 green 0.6 blue 0.6 
//light_source { <3*LL, 0*LL, 2*LL> color red 0.6 green 0.6 blue 0.6 
    area_light <0, 1, 0>, <0, 0, 1>, 3, 3 // area light along yz-plane with size 3x3 lamps
    adaptive 1
    jitter
}


#declare T_BoxSurf =
  texture{pigment{color White transmit 0.9 } finish{phong .8 phong_size 200 reflection 0.02}} //use for black background
// texture{pigment{color Col_Glass_Clear} finish{phong .8 phong_size 200 reflection 0.02} } //use for white background

// Box with some thickness, describing the cell
#declare Box =
difference{
box{ //outer box
<-0.5*LL2,-0.5*LL2,-0.5*LL2>, <0.5*LL2,0.5*LL2,0.5*LL2>
}
box{ //inner box
  <-0.5*LL,-0.5*LL,-0.5*LL>, <0.5*LL,0.5*LL,0.5*LL>
}
texture{ T_BoxSurf }
interior {ior iorCylinder} //to make the inside of the cylinder glass like 
//interior { I_Glass1 } //use for white background
} //end of difference

//Make a box that describes the fluid
#declare Box2 =
box{
  <-0.5*LL,-0.5*LL,-0.5*LL>, <0.5*LL,0.5*LL,0.5*LL>
  texture{ T_BoxSurf }
// DO NOT USE INTERIOR BECAUSE it "multiplies" the visible vortices 
//  interior {ior iorCylinder} //to make the inside of the cylinder glass like 
  hollow on
}

object{
//Box  //Thick walls
Box2 //Thin walls
}

#include "mesh.pov" // Include the vortex data. Note the vortex object should have a name "Vortices"

object {Vortices}


object{
ColorBar
scale <10,10,10>
//rotate -90*x
translate 7*y
translate -4*x
}


