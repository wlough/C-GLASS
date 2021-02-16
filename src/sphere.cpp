#include "cglass/sphere.hpp"

Sphere::Sphere(unsigned long seed) : Object(seed) {
  shape_ = shape::sphere;
}