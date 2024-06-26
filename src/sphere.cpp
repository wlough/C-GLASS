#include "cglass/sphere.hpp"

Sphere::Sphere(unsigned long seed) : Object(seed) {
  shape_ = shape::sphere;
}

//Set the neighboring receptors
void Sphere::SetNeighbors(Sphere* prev_r, Sphere* next_r) {
  prev_r_ = prev_r;
  next_r_ = next_r;
}

//Set point cover object sphere is on
void Sphere::SetPCObjectForSphere(Object* pc_object) {
  pc_object_ = pc_object;
}

//Get orientation of point cover object
Object* Sphere::GetPCObjectForSphere() {
  return pc_object_;
}

//Set how long the sphere is along the filament
void Sphere::SetSphereS(double s) {
  s_ = s;
}

//Get how far the sphere is along the filament
double Sphere::GetSphereS() {
  return s_;
}

//Set the spacing between sites
void Sphere::SetStepSize(double step_size) {
  step_size_ = step_size;
} 

//Get the spacing between sites 
double Sphere::GetStepSize(){
  return step_size_;
}

//Get the plus neighbor of current receptor
Sphere* Sphere::GetPlusNeighbor() {
  return next_r_;
}

//Get the minus neighbor of current receptor
Sphere* Sphere::GetMinusNeighbor() {
  return prev_r_;
}
