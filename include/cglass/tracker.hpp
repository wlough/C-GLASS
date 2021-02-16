#ifndef _CGLASS_TRACKER_H_
#define _CGLASS_TRACKER_H_

class Tracker {
  private:
    #ifdef TRACK
      // Unbound-to-Single variables
      int n_us_ = 0; // number of crosslinks that went from unbound to single
      double n_us_predicted_ = 0.0; // expected crosslinks that go from unbound to single

      // Single-to-Unbound variables
      int n_su_ = 0; // number of crosslinks that went from single to unbound
      double n_su_predicted_ = 0.0; // expected crosslinks that go from single to unbound

      // Single-to-Double variables
      int n_sd_ = 0; // number of crosslinks that went from single to double
      double n_sd_predicted_ = 0.0; // expected crosslinks that go from single to double

      // Double-to-Single variables
      int n_ds_ = 0; // number of crosslinks that went from double to single
      double n_ds_predicted_ = 0.0; // expected crosslinks that go from double to single
    #endif
  public:
    // consider a S->D crosslink transition
    void TrackUS(double n_expected) {
      #ifdef TRACK
        n_us_predicted_ += n_expected;
      #endif
    }
    void TrackSU(double p_su) {
      #ifdef TRACK
        n_su_predicted_ += p_su;
      #endif
    }
    void TrackSD(double p_sd) {
      #ifdef TRACK
        n_sd_predicted_ += p_sd;
      #endif
    }
    void TrackDS(double p_ds) {
      #ifdef TRACK
        n_ds_predicted_ += p_ds;
      #endif
    }
    // track binding from S->D
    void BindUS(int n) { 
      #ifdef TRACK
        n_us_+=n; 
      #endif
    }
    void UnbindSU() { 
      #ifdef TRACK
        n_su_++;
      #endif
    }
    void BindSD() { 
      #ifdef TRACK
        n_sd_++;
      #endif
    }
    void UnbindDS() {  
      #ifdef TRACK
        n_ds_++;
      #endif
    }
    void PrintTrack() {
      #ifdef TRACK
        Logger::Track("n_us predicted: %f", n_us_predicted_);
        Logger::Track("n_us: %d", n_us_);
        Logger::Track("n_su predicted: %f", n_su_predicted_);
        Logger::Track("n_su: %d", n_su_);
        Logger::Track("n_sd predicted: %f", n_sd_predicted_);
        Logger::Track("n_sd: %d", n_sd_);
        Logger::Track("n_ds predicted: %f", n_ds_predicted_);
        Logger::Track("n_ds: %d", n_ds_);
      #endif
    }
};

#endif // _CGLASS_TRACKER_H_