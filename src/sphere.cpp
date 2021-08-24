#include "cglass/sphere.hpp"

Sphere::Sphere(unsigned long seed) : Object(seed) {
  shape_ = shape::sphere;
}
void Sphere::SetSphereReceptor(Object* sphere_receptor){
     sphere_receptor_=sphere_receptor;
}	
Object* Sphere::GetSphereReceptor(){
return sphere_receptor_;
}
void Sphere::SetNeighbors(Object* prevN, Object* nextN) {
  prevN_ = prevN;
  nextN_ = nextN;

}
Object* Sphere::GetPlusNeighbor() {
	    //Logger::Warning("nextN_ID %i", nextN_->GetOID());
	        return nextN_;
}
Object* Sphere::GetMinusNeighbor() {
	    //Logger::Warning("nextN_ID %i", nextN_->GetOID());
	        return prevN_;
}


