#include "cglass/point_cover.hpp"
#include "cglass/receptor.hpp"

PointCover::PointCover(unsigned long seed) : Object(seed) {
  comp_type_ = comp_type::point_cover;
}

void PointCover::AddSpherePtr(Receptor *s) {
  sphere_ptrs_.push_back(s);
  sphere_ptrs_.back()->SetColor(color_, draw_);
  sphere_ptrs_.back()->SetCompID(GetCompID());
}