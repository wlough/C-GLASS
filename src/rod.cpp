#include "cglass/rod.hpp"

Rod::Rod(unsigned long seed) : Object(seed) {
  shape_ = shape::rod;
}