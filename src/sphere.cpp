#include "cglass/sphere.hpp"

Sphere::Sphere(unsigned long seed) : Object(seed) {
  shape_ = shape::sphere;
}

//Set the neighboring receptors
void Sphere::SetNeighbors(Object* prev_r, Object* next_r) {
  prev_r_ = prev_r;
  next_r_ = next_r;
}

//Set point cover object sphere is on
void Sphere::SetPCObjectForSphere(Object* pc_object) {
  pc_object_ = pc_object;
  Logger::Warning("pc_object for sphere set to %i", pc_object_->GetOID());
}

//Get orientation of point cover object
Object* Sphere::GetPCObjectForSphere() {
  return pc_object_;
}

//Get the plus neighbor of current receptor
Object* Sphere::GetPlusNeighbor() {
 return next_r_;
}

//Get the minus neighbor of current receptor
Object* Sphere::GetMinusNeighbor() {
 return prev_r_;
}
