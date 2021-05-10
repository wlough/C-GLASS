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
                                 " n_steps: 10000,"
                                 " delta: 0.001,"
                                 " n_dim: 3,"
                                 " n_periodic: 3,"
                                 " system_radius: 20,"
                                 " br_bead: {num: 100}}";

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
    SECTION("Test diffusion relations") {
      REQUIRE(TestInit(run_opts));
      REQUIRE(TestSim());
    }
  }
};

TEST_CASE("Test simulation manager") {
  UnitTest<SimulationManager> test;
  test.TestSimManager();
  test.TestDiffusion();
}
