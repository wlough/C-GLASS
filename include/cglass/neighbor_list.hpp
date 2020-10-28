#include "auxiliary.hpp"
#include "object.hpp"
#include "bond.hpp"
#include "site.hpp"
#include <mutex>

/* A data structure that is used to hold a list of particles that are nearby the
 * owner of the list in a thread-safe way. */
class NeighborList {
private:
  std::mutex mtx_;
  std::vector<Object *> nlist_;
  std::vector<Bond *> nlist_bond_;
  std::vector<Site *> nlist_site_;

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
    switch(obj->GetType()) { 
      case obj_type::bond:
        nlist_bond_.push_back(dynamic_cast<Bond*>(obj));
        break;
      case obj_type::site:
        nlist_site_.push_back(dynamic_cast<Site*>(obj));
        break;
      default:
        break;
    }
  }
  const Object *const *GetNeighborListMem() { return &nlist_[0]; }
  
  const std::vector<Bond*>& GetNeighborListMemBonds() {
    return nlist_bond_;
  }

  const std::vector<Site*>& GetNeighborListMemSites() {
    return nlist_site_;
  }
  void Clear() { 
    nlist_.clear();
    nlist_bond_.clear();
    nlist_site_.clear();
   }
  const int NNeighbors() const { return nlist_.size(); }

  const int NNeighborsSite() const { return nlist_site_.size(); }

  const int NNeighborsBond() const { return nlist_bond_.size(); }

  Object *GetNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_.size()) {
      Logger::Error("Invalid index received in class NeighborList");
    }
    return nlist_[i_neighbor];
  }

  Site *GetSiteNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_site_.size()) {
      Logger::Error("Invalid index received in class NeighborList");
    }
    return nlist_site_[i_neighbor];
  }

  Bond *GetBondNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_bond_.size()) {
      Logger::Error("Invalid index received in class NeighborList");
    }
    return nlist_bond_[i_neighbor];
  }
};
