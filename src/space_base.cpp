#include "cglass/space_base.hpp"

double SpaceBase::BoundaryArea() const{
  switch (type) {
    case +boundary_type::none: {
      Logger::Error("Boundary area calculation requires a boundary.");
      break;
    }
    case +boundary_type::sphere: {
      if (n_dim == 2) {
        return 2.0 * M_PI * radius;
      } else {
        return 4 * M_PI * SQR(radius);
      }
    }
    case +boundary_type::protrusion: {
      if (n_dim == 2) {
        Logger::Error("Protrusion not set up for 2d");
      } else {
        return 4* M_PI * SQR(radius) + M_PI * SQR(pro_radius) * pro_length;
      }
    }
    case +boundary_type::box: {
      if (n_dim == 2) {
        return 8.0 * radius;
      } else {
        return 24.0 * SQR(radius);
      }
    }
    case +boundary_type::budding: {
      double R = radius;
      double r = bud_radius;
      double d = bud_height;
      if (n_dim == 2) {
        // arc length is 2*r*(pi-theta) where theta is intersect angle
        return 2.0 * (r * (M_PI - acos((SQR(d) + SQR(r) - SQR(R)) / (2.0 * d * r))) 
               + R * (M_PI - acos((SQR(d) - SQR(r) + SQR(R)) / (2.0 * d * R))));
      } else {
        // segment area is 2*pi*r^2*(1+cos(theta)) where theta is intersect angle
        return 2.0 * M_PI * (SQR(r) * (1 + ((SQR(d) + SQR(r) - SQR(R)) / (2.0 * d * r))) 
               + SQR(R) * (1 + ((SQR(d) - SQR(r) + SQR(R)) / (2.0 * d * R))));
      }
    }
    case +boundary_type::wall: {
      Logger::Error("Cannot calculate area for wall boundary.");
      break;
    }
    default: 
      Logger::Error("Boundary type not recognized in Space Base.");
  }
  return -1;
}
