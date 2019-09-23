  #include "colors.inc"
  #fopen MyFile "test.log" read
  #while (defined(MyFile))
   #declare Var2=<0,0,0>;
     #declare Var1=<0,0,0>;
     #read (MyFile,Var1,Var2)
  cylinder {
    Var1,     // Center of one end
    Var2,     // Center of other end
    .0004            // Radius
   pigment { Red }
   finish { ambient 0.2 diffuse 0.9 phong 1 }
  }
  #end
 background {Black}
//box
//{
//< -0.05, -0.05, -0.05> < 0.05, 0.05, 0.05>
//pigment {color rgbf < 0.8, .9, 0.9,.99>}
//   pigment { White filter .9 }
//finish {
//      ambient 0
//      diffuse 0.
//      reflection .1
//      specular 1
//      roughness .001
//   }
//}
//====== Cube formed by tubes =====
#declare R = 0.0005; //radius of tubes
#declare BigCube1 =
union{
// 8 Corners
sphere{<-0.05,-0.05,-0.05>,R}
sphere{< 0.05,-0.05,-0.05>,R}
sphere{<-0.05,-0.05, 0.05>,R}
sphere{< 0.05,-0.05, 0.05>,R}
sphere{<-0.05, 0.05,-0.05>,R}
sphere{< 0.05, 0.05,-0.05>,R}
sphere{<-0.05, 0.05, 0.05>,R}
sphere{< 0.05, 0.05, 0.05>,R}
// 4 in x direction
cylinder {<-0.05,-0.05,-0.05>,< 0.05,-0.05,-0.05>,R}
cylinder {<-0.05,-0.05, 0.05>,< 0.05,-0.05, 0.05>,R}
cylinder {<-0.05, 0.05,-0.05>,< 0.05, 0.05,-0.05>,R}
cylinder {<-0.05, 0.05, 0.05>,< 0.05, 0.05, 0.05>,R}
// 4 in y direction
cylinder {<-0.05,-0.05,-0.05>,<-0.05, 0.05,-0.05>,R}
cylinder {<-0.05,-0.05, 0.05>,<-0.05, 0.05, 0.05>,R}
cylinder {< 0.05,-0.05,-0.05>,< 0.05, 0.05,-0.05>,R}
cylinder {< 0.05,-0.05, 0.05>,< 0.05, 0.05, 0.05>,R}
// 4 in z direction
cylinder {<-0.05,-0.05,-0.05>,<-0.05,-0.05, 0.05>,R}
cylinder {<-0.05, 0.05,-0.05>,<-0.05, 0.05, 0.05>,R}
cylinder {< 0.05,-0.05,-0.05>,< 0.05,-0.05, 0.05>,R}
cylinder {< 0.05, 0.05,-0.05>,< 0.05, 0.05, 0.05>,R}
texture{pigment{color rgb<.9,.9,.9>}
        finish{ diffuse 0.9 phong 1}}
}//-- End of wireframed cube -------

//------------- Draw it  -----------
object{BigCube1}
//------------------------------ end
camera {
 sky<0,0,1>
  location <0.2, -0.2, 0.1>
  look_at <0, 0, 0>
  angle 43
}
//light_source { <0.0, -0.02, -0.2> color White}
//light_source { <-0.2, 0.2, 0.2> color White}
//light_source { <-0.5, 0.5, 0.5> color White}
  light_source {
    <1, -1, -1>
    color White
    area_light <1, 0, 0>, <0, 0, 1>, 5, 5
    adaptive 1
    jitter
  }

