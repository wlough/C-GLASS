
#include <iostream>
#include <stdlib.h>
#include "parse_flags.h"
#include "simulation_manager.h"

/*************************
   ::SimCORE Main::
   Parse commandline flags and start simulation manager
**************************/
int main(int argc, char *argv[]) {

  // Parse input flags, see parse_flags.h for documentation
  run_options run_opts = parse_opts(argc, argv);
  if (!run_opts.f_flag && !run_opts.debug) {
    std::cout << "  ERROR: No parameter file given!\n";
    show_help_info(argv[0]);
    exit(1);
  }
  SimulationManager sim;

  if (run_opts.debug) {
    // Debug mode should be very hands-off
    sim.DebugMode();
  }

  if (run_opts.test) {
    // Run the test mode operation
    debug_trace = true;
  } 

  else {
    // Initialize param_file, rng, run_name, n_runs (if in param_file)
    sim.InitManager(run_opts.param_file);
    // Prefer command-line options over param values for n_runs, run_name
    if (run_opts.n_flag)
      sim.SetNRuns(run_opts.n_runs);
    if (run_opts.r_flag)
      sim.SetRunName(run_opts.run_name);
    if (run_opts.m_flag)
      // Recreate Movie of past simulations using posit file
      sim.RunMovieManager(run_opts.posit_files);
    else if (run_opts.a_flag)
      // Analyze existing posit files
      sim.RunAnalyses(run_opts.posit_files);
    else
      // Main control function
      sim.RunManager();
  }

  return 0;
}

