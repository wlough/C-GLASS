#ifndef _CGLASS_TRACKER_H_
#define _CGLASS_TRACKER_H_

class Tracker {
  private:
    int n_sd_; // number of crosslinks that went from single to double
    int n_sd_considered_; // number of crosslinks that considered going from single to double
    int n_sd_predicted_; // expected crosslinks that go from single to double
  public:
    // consider a S->D crosslink transition
    void TrackSD(double p_sd) {
      n_sd_considered_++;
      n_sd_predicted_ += p_sd;
    }
    // track binding from S->D
    void BindSD() { n_sd_++; }
    void PrintTrack() {
      std::cout << "n_sd predicted: " << n_sd_predicted_ << std::endl;
      std::cout << "n_sd: " << n_sd_ << std::endl;
    }
};

#endif // _CGLASS_TRACKER_H_