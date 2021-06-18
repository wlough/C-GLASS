// Enable unit testing always
#ifndef TESTS
#define TESTS
#endif

#include <cglass/simulation_manager.hpp>
#include "catch.hpp"

template <> class UnitTest<SimulationManager> {
private:
  SimulationManager mgr_;
  const int n_xlink_sd_ = 2;
  const int n_recep_sd_ = 3;
  std::string multi_params_ = "{run_name: test,"
                              " seed: 314159,"
                              " n_runs: 3,"
                              " n_random: 2,"
                              " n_steps: 42,"
                              " n_dim: [V, 2, 3],"
                              " delta: [R, 0.1, 0.001],"
                              " system_radius: 1}";
  std::string diffusion_params_ = "{run_name: diffusion_test,"
                                  " seed: 314159,"
                                  " n_steps: 1000,"
                                  " n_periodic: 3,"
                                  " system_radius: 20,"
                                  " potential: none,"
                                  " br_bead: {num: 100,spec_flag: true,n_spec: 500}}";
  std::string SD_params_ = "{run_name: sd_test,"
                           " seed: 314159,"
                           " n_steps: 100000,"
                           " n_dim: 2,"
                           " no_midstep: true,"
                           " system_radius: 20,"
                           " rigid_filament: [{name: rig1, num: 1, length: 5, stationary_flag: true,"
                             " insertion_type: custom, insert_file: fil1.yaml},"
                             " {name: rig2, num: 1, length: 5, stationary_flag: true,"
                             " insertion_type: custom, insert_file: fil2.yaml}],"
                           " crosslink: {concentration: 1, begin_with_bound_crosslinks: 2, k_spring: 2,"
                             " anchors: [bind_file: rig_bind.yaml]},"
                           " receptor: [{name: rig1_receptor, num: 2, diameter: 0.9, component: rig1, insertion_type: grid},"
                             " {name: rig2_receptor, num: 3, diameter: 0.9, component: rig2, spec_flag: true, n_spec: 20,"
                             " insertion_type: grid}]}";
  std::string fil1_insert_ = "rigid_filament: [[[0,0,0],[0.99619469809,0.08715574274,0]]]";
  std::string fil2_insert_ = "rigid_filament: [[[0,1.5,0],[1,0,0]]]";
  std::string rig_bind_ = "{rig1_receptor: {k_on_s: 5, k_off_s: 0, k_on_d: 0, k_off_d: 0,"
                          " density_type: surface, bind_site_density: 0.4, single_occupancy: true},"
                          " rig2_receptor: {k_on_s: 0, k_off_s: 0, k_on_d: 100, k_off_d: 10,"
                          " density_type: surface, bind_site_density: 0.4, single_occupancy: true}}";
  // Helper function to process time header in txt file, returns time
  double ProcessTimeHeader(std::ifstream& filein) {
    std::string line;
    double t = -1;
    getline(filein, line);
    t = stod(line.substr(7, line.size()-7));
    for (int i = 0; i < 2; ++i) getline(filein, line);
    return t;
  }
  // Helper function to determine theoretical % bound from partition function
  void CalcTheoreticalFraction(double* fraction) {
    // Parameters used in simulation
    double kon = 100.0;
    double koff = 10.0;
    double k = 2.0;
    double density = 0.4;
    double radius = 0.45;
    double small = 1e-4;
    // Derived parameters
    double kon_true, bind_energy;
    double energy[n_xlink_sd_][n_recep_sd_];
    double Z = 1.0; // partition function
    // Calculate bind volume as KMC does
    double bind_vol = 4 * M_PI * boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
                      [&](double s){return SQR(s) * exp(-k*SQR(s)/2);}, 0, sqrt(-2*log(small)/k), 
                      10, 1e-6, nullptr);
    // Calculate true on-rate using binding volume, density, site size
    kon_true = kon * density * 4 * M_PI * SQR(radius) / bind_vol;
    bind_energy = log(kon_true/koff);
    // Fill a vector of energies for each interaction
    for (int i = 0; i < n_xlink_sd_; ++i) {
      for (int j = 0; j < n_recep_sd_; ++j) {
        // position along rigid_filament
        double xlink_s = (-2.5 + 2.5*(i+0.5));
        // r_sqr is (xlink_pos_i_x-site_pos_j_x)^2+(xlink_pos_i_y-site_pos_j_y)^2
        double r_sqr = SQR(xlink_s*cos(5*M_PI/180)-(-2.5 + 5.0/3*(j+0.5)))
                     + SQR(xlink_s*sin(5*M_PI/180)-1.5);
        energy[i][j] = k*r_sqr/2 - bind_energy;
      }
    }
    // Calculate the partition function & theoratical fractions 
    // (the receptors can only bind to one xlink)
    for (int i = 0; i < n_xlink_sd_; ++i) {
      for (int j = 0; j < n_recep_sd_; ++j) {
        double boltzmann = exp(-energy[i][j]);
        Z += boltzmann;
        fraction[j] += boltzmann;
        // Start at i+1 to eliminate double-counts
        for (int ib = i + 1; ib < n_xlink_sd_; ++ib) {
          for (int jb = 0; jb < n_recep_sd_; ++jb) {
            if (ib != i && jb != j) {
              boltzmann = exp(-energy[i][j]-energy[ib][jb]);
              Z += boltzmann;
              fraction[j] += boltzmann;
              fraction[jb] += boltzmann;
            }
          }
        }
      }
    }
    for (int j = 0; j < n_recep_sd_; ++j) {
      fraction[j] /= Z;
    }
  } // Helper function to determine fraction bound
  void CalcSimulationFraction(double* fraction) {
    // Parameters used in simulation
    std::ifstream fin("sd_test_receptor_rig2_receptorSpec.txt", std::ios::in);
    std::string line;
    int total_points = 0;
    // Ignore header
    for (int i = 0; i < 2; ++i) getline(fin, line);
    for (int i_step = 0; i_step < 5000; ++i_step) {
      ProcessTimeHeader(fin);
      for (int rec_index = 0; rec_index < n_recep_sd_; ++rec_index) {
        if (!getline(fin, line)) {
          std::cerr << "Output file incomplete in CalcSimulationFraction" << std::endl;
          fraction[rec_index] = -1;
          return;
        }
        // Allow to reach equilibrium before calculation
        if (i_step >= 100) {
          total_points++;
          // Add number of bound xlinks
          fraction[rec_index] += stod(line);
        }
      }
    }
    for (int rec_index = 0; rec_index < n_recep_sd_; ++rec_index) {
      fraction[rec_index] /= (total_points / n_recep_sd_);
    }
    fin.close();
  }
  // Helper function for determining diffusion from spec file
  double CalcDiffusionCoef() {
    std::ifstream fin("diffusion_test_br_bead_speciesSpec.txt", std::ios::in);
    std::stringstream sin;
    std::string line;
    int nbeads = 100;
    std::vector<std::vector<double> > r0(nbeads,std::vector<double>(3));
    double rij;
    double dr2 = 0; // squared displacement
    double ti, tf;
    // Process headers
    for (int i = 0; i < 2; ++i) getline(fin, line);
    ti = ProcessTimeHeader(fin);
    // Get initial positions of beads
    for (int i = 0; i < nbeads; ++i) {
      if (!getline(fin, line)) return -1.0;
      sin = std::stringstream(line);
      for (int j = 0; j < 3; ++j) {
        sin >> r0[i][j];
      }
    }
    tf = ProcessTimeHeader(fin);
    // Calculate msd
    for (int i = 0; i < nbeads; ++i) {
      if (!getline(fin, line)) return -1.0;
      sin = std::stringstream(line);
      for (int j = 0; j < 3; ++j) {
        sin >> rij;
        dr2 += (rij-r0[i][j])*(rij-r0[i][j]);
      }
    }
    fin.close();
    // Alpha = 3 for 3D, D should be 1 in sim units
    return dr2/(3*nbeads*(tf-ti)*2);
  }
public:
  void TestSimManager() {
    std::string fname = "test_params.yaml";
    std::ofstream(fname, std::ios::out) << multi_params_;
    run_options run_opts;
    run_opts.param_file = fname;
    SECTION("Simulation manager init and run simulations") {
      REQUIRE(TestInit(run_opts));
      TestParams();
      REQUIRE(TestSim());
    }
  }
  bool TestInit(run_options run_opts) {
    try {
      mgr_.InitManager(run_opts);
    } catch (int e) {
      printf("Exception occurred with No. %d", e);
      return false;
    }
    return true;
  }
  void TestParams() {
    REQUIRE(mgr_.run_name_.compare("test") == 0);
    REQUIRE(mgr_.n_runs_ == 3);
    REQUIRE(mgr_.n_random_ == 2);
    // We have two parameters for n_dim and two random variables, 2*2=4
    REQUIRE(mgr_.n_var_ == 4);
  }
  bool TestSim() {
    try {
      mgr_.RunManager();
    } catch (int e) {
      printf("Exception occurred with No. %d", e);
      return false;
    }
    return true;
  }
  void TestDiffusion() {
    std::string fname = "diffusion_test_params.yaml";
    std::ofstream(fname, std::ios::out) << diffusion_params_;
    run_options run_opts;
    run_opts.param_file = fname;
    SECTION("Run diffusion test") {
      REQUIRE(TestInit(run_opts));
      REQUIRE(TestSim());
    }
    run_opts.convert = true;
    // Reset to clear cell lists, etc.
    SECTION("Convert diffusion test results") {
      REQUIRE(TestInit(run_opts));
      REQUIRE(TestSim());
    }
    SECTION("Analyze diffusion test results") {
      // Use 3 sigma to test D (should be 1 in sim units)
      REQUIRE(CalcDiffusionCoef()<1.2);
      REQUIRE(CalcDiffusionCoef()>0.8);
    }
  }
  void TestSD() {
    std::string fname = "sd_test_params.yaml";
    std::ofstream(fname, std::ios::out) << SD_params_;
    std::ofstream("fil1.yaml", std::ios::out) << fil1_insert_;
    std::ofstream("fil2.yaml", std::ios::out) << fil2_insert_;
    std::ofstream("rig_bind.yaml", std::ios::out) << rig_bind_;
    run_options run_opts;
    run_opts.param_file = fname;
    SECTION("Run single-double (SD) binding test") {
      REQUIRE(TestInit(run_opts));
      REQUIRE(TestSim());
    }
    run_opts.convert = true;
    // Reset to clear cell lists, etc.
    SECTION("Convert SD test results") {
      REQUIRE(TestInit(run_opts));
      REQUIRE(TestSim());
    }
    SECTION("Analyze SD test results") {
      double th_fraction[3] = {0,0,0};
      double sim_fraction[3] = {0,0,0};
      // 3 sigma determined by testing
      double three_sigma[3] = {0.028,0.031,0.056};
      CalcTheoreticalFraction(th_fraction);
      CalcSimulationFraction(sim_fraction);
      for (int rec_index = 0; rec_index < n_recep_sd_; ++ rec_index) {
        REQUIRE(sim_fraction[rec_index] > (th_fraction[rec_index]-three_sigma[rec_index]));
        REQUIRE(sim_fraction[rec_index] < (th_fraction[rec_index]+three_sigma[rec_index]));
      }
    }
  }
};

TEST_CASE("Test simulation manager") {
  UnitTest<SimulationManager> test;
  test.TestSimManager();
  test.TestDiffusion();
  test.TestSD();
}
