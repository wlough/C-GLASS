#ifndef _CGLASS_SPHERE_H_
#define _CGLASS_SPHERE_H_

#include "object.hpp"

// General class of spherically shaped objects; all spherical objects can only be bound once 
// to an xlink (in their centers). KMC is constructed with TSphere=Sphere in C-GLASS.
class Sphere : public Object {
    protected:
      bool induces_catastrophe_ = false; // Whether the object induces catastrophe upon double binding
    private:  
      Sphere* prev_r_ = nullptr; // The previous receptor (in - direction) on the filament 
      Sphere* next_r_ = nullptr; // The next receptor (in + direction) on the filament
      Object* pc_object_ = nullptr; // point cover object sphere is on
      double step_size_ = 0;
    public:
      Sphere(unsigned long seed);
      bool InducesCatastrophe() { return induces_catastrophe_; }
      void SetNeighbors(Sphere* prevN, Sphere* nextN); //Set neighboring receptors on filament
      void SetPCObjectForSphere(Object* pc_object);
      void SetStepSize(double step_size);
      double GetStepSize();
      Object* GetPCObjectForSphere();
      Sphere* GetPlusNeighbor();
      Sphere* GetMinusNeighbor();
};
#endif // _CGLASS_SPHERE_H_
