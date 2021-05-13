// Enable unit testing always
#ifndef TESTS
#define TESTS
#endif

#include <cglass/simulation_manager.hpp>
#include "catch.hpp"

template <> class UnitTest<SimulationManager> {
private:
  SimulationManager mgr_;
  std::string multi_params_ = "{run_name: test,"
                              " seed: 314159,"
                              " n_runs: 3,"
                              " n_random: 2,"
                              " n_steps: 42,"
                              " n_dim: [V, 2, 3],"
                              " delta: [R, 0.1, 0.001],"
                              " system_radius: 1}";
  std::string diffusion_params_ = "{run_name: test,"
                                 " seed: 314159,"
                                 " n_steps: 1000,"
                                 " delta: 0.001,"
                                 " n_dim: 3,"
                                 " n_periodic: 3,"
                                 " system_radius: 20,"
                                 " potential: none,"
                                 " br_bead: {num: 100,spec_flag: true,n_spec: 500,midstep: false}}";
  // Helper function to process time header in txt file, returns time
  double ProcessTimeHeader(std::ifstream& filein) {
    std::string line;
    double t = -1;
    getline(filein, line);
    t = stod(line.substr(7, line.size()-7));
    for (int i = 0; i < 2; ++i) getline(filein, line);
    return t;
  }
  // Helper function for determining diffusion from spec file
  double CalcDiffusionCoef() {
    std::ifstream fin("test_br_bead_speciesSpec.txt", std::ios::in);
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
    std::string fname = "test_params.yaml";
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
};

TEST_CASE("Test simulation manager") {
  UnitTest<SimulationManager> test;
  test.TestSimManager();
  test.TestDiffusion();
}
