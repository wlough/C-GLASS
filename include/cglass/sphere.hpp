#ifndef _CGLASS_SPHERE_H_
#define _CGLASS_SPHERE_H_

#include "object.hpp"
//#include "receptor.hpp"
// General class of spherically shaped objects; all spherical objects can only be bound once 
// to an xlink (in their centers). KMC is constructed with TSphere=Sphere in C-GLASS.
class Sphere : public Object {
    protected:
      bool induces_catastrophe_ = false; // Whether the object induces catastrophe upon double binding
      Object* sphere_receptor_ =nullptr;
      Object* prevN_ = nullptr;
      Object* nextN_ = nullptr;
    public:
      Sphere(unsigned long seed);
      bool InducesCatastrophe() { return induces_catastrophe_; }
      void SetSphereReceptor(Object* sphere_receptor);
      void SetNeighbors(Object* prevN, Object* nextN);
      Object* GetPlusNeighbor();
      Object* GetMinusNeighbor();
      Object* GetSphereReceptor();
      
};
#endif // _CGLASS_SPHERE_H_
