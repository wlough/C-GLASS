#ifndef _CGLASS_TRACKER_H_
#define _CGLASS_TRACKER_H_

class Tracker {
  private:
    // Single-to-Double variables
    int n_sd_ = 0; // number of crosslinks that went from single to double
    int n_sd_considered_ = 0; // number of crosslinks that considered going from single to double
    double n_sd_predicted_ = 0.0; // expected crosslinks that go from single to double

    // Single-to-Unbound variables
    int n_su_ = 0; // number of crosslinks that went from single to unbound
    int n_su_considered_ = 0; // number of crosslinks that considered going from single to unbound
    double n_su_predicted_ = 0.0; // expected crosslinks that go from single to unbound
  public:
    // consider a S->D crosslink transition
    void TrackSD(double p_sd) {
      n_sd_considered_++;
      n_sd_predicted_ += p_sd;
    }
    void TrackSU(double p_su) {
      n_su_considered_++;
      n_su_predicted_ += p_su;
    }
    // track binding from S->D
    void BindSD() { n_sd_++; }
    void UnbindSU() { n_su_++; }
    void PrintTrack() {
      std::cout << "n_sd predicted: " << n_sd_predicted_ << std::endl;
      std::cout << "n_sd: " << n_sd_ << std::endl;
      std::cout << "n_su predicted: " << n_su_predicted_ << std::endl;
      std::cout << "n_su: " << n_su_ << std::endl;
    }
};

#endif // _CGLASS_TRACKER_H_