
/* A data structure that is used to hold a list of particles that are nearby the
 * owner of the list in a thread-safe way. */
#include "auxiliary.hpp"
#include "bond.hpp"
#include "receptor.hpp"
#include <mutex>

class NeighborList {
private:
  std::mutex mtx_;
  std::vector<Object *> nlist_;
  std::vector<const Bond *> nlist_rod_;
  std::vector<const Receptor *> nlist_sphere_;

public:
  NeighborList() {}
  ~NeighborList() { Clear(); }
  NeighborList(const NeighborList &that) { this->nlist_ = that.nlist_; }
  NeighborList &operator=(NeighborList const &that) {
    this->nlist_ = that.nlist_;
    return *this;
  }
  void AddNeighbor(Object *obj) {
    std::lock_guard<std::mutex> lk(mtx_);
    nlist_.push_back(obj);
    switch (obj->GetShape()) {
    case shape::rod:
      nlist_rod_.push_back(dynamic_cast<const Bond *>(obj));
      break;
    case shape::sphere:
      if (obj->GetNAnchored() == 0) {
        nlist_sphere_.push_back(dynamic_cast<const Receptor *>(obj));
      }
      break;
    default:
      break;
    }
  }
  const Object *const *GetNeighborListMem() { return &nlist_[0]; }

  const std::vector<const Bond *> &GetNeighborListMemRods() {
    return nlist_rod_;
  }

  const std::vector<const Receptor *> &GetNeighborListMemSpheres() {
    return nlist_sphere_;
  }
  void Clear() {
    nlist_.clear();
    nlist_rod_.clear();
    nlist_sphere_.clear();
  }
  const int NNeighbors() const { return nlist_.size(); }

  const int NNeighborsSphere() const { return nlist_sphere_.size(); }

  const int NNeighborsRod() const { return nlist_rod_.size(); }

  Object *GetNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_.size()) {
      Logger::Error("Invalid index received in class NeighborList");
    }
    return nlist_[i_neighbor];
  }

  Receptor *GetSphereNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_sphere_.size()) {
      Logger::Error("Invalid index received in class NeighborList");
    }
    return const_cast<Receptor *>(nlist_sphere_[i_neighbor]);
  }

  Bond *GetRodNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_rod_.size()) {
      Logger::Error("Invalid index received in class NeighborList");
    }
    return const_cast<Bond *>(nlist_rod_[i_neighbor]);
  }
};
