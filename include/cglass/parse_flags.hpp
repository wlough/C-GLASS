#ifndef _CGLASS_PARSE_FLAGS_H_
#define _CGLASS_PARSE_FLAGS_H_

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include "version_generator.hpp"

// Run parameters that needs to be carried forward should be stored here.
struct run_options {
  int n_runs = 1;
  bool debug = false;
  bool run_name_flag = false;
  bool n_run_flag = false;
  bool graphics_flag = false;
  bool make_movie = false;
  bool analysis_flag = false;
  bool use_posits = false;
  bool load_checkpoint = false;
  bool reduce_flag = false;
  int n_graph = 100;
  int reduce_factor = 1;
  bool blank_flag = false;
  bool auto_graph = false;
  bool with_reloads = false;
  bool single_frame = false;
  bool convert = false;
  std::string param_file;
  std::string run_name = "sc";
};

/* NOTE TO THE FUTURE: Need to update n_flags, PARSE_OPTS function
   case-switches, as well as any flag descriptions in the static
   string array below when adding new flags. */

// Define flags here
static const int n_flags = 16;
static struct option long_options[] = {{"help", no_argument, 0, 'h'},
                                       {"version", no_argument, 0, 'v'},
                                       {"debug", no_argument, 0, 'd'},
                                       {"run-name", required_argument, 0, 'r'},
                                       {"n-runs", required_argument, 0, 'n'},
                                       {"graphics", required_argument, 0, 'g'},
                                       {"movie", no_argument, 0, 'm'},
                                       {"auto-graph", no_argument, 0, 'G'},
                                       {"analysis", no_argument, 0, 'a'},
                                       {"posit", no_argument, 0, 'p'},
                                       {"load", no_argument, 0, 'l'},
                                       {"reduce", required_argument, 0, 'R'},
                                       {"blank", no_argument, 0, 'b'},
                                       {"with-reloads", no_argument, 0, 'w'},
                                       {"single-frame", no_argument, 0, 'M'},
                                       {"convert", no_argument, 0, 'c'},
                                       {0, 0, 0, 0}};

// Descriptions for flags
static const std::string desc[n_flags][2] = {
    {"show this menu", "none"},
    {"show C-GLASS version number", "none"},
    {"run program in debug mode", "none"},
    {"where rname is the name of a run session (for organizing batch jobs)",
     "rname"},
    {"where num is the number of independent runs to perform with the given "
     "parameters",
     "num"},
    {"force to run with graphics from posit/spec files, with n_graph = ngraph",
     "ngraph"},
    {"run simulation to make movie using spec files with same run-name",
     "none"},
    {"if running with graphics, do not wait for human input at graphics "
     "initialization",
     "none"},
    {"run analysis on spec files with same run-name for parameter species",
     "none"},
    {"use posit files for movies/analysis rather than spec files", "none"},
    {"run simulation starting from checkpoint files with same run-name for "
     "parameter species",
     "none"},
    {"reduce output file resolution by a factor of reduce_factor",
     "reduce_factor"},
    {"write all necessary parameter files without running any simulations",
     "none"},
    {"run analysis using any existing spec files from reloaded runs", "none"},
    {"runs movie, but only generates a single bitmap image to record the final "
     "system state",
     "none"},
    {"converts posit and spec files to text files",
     "none"}};

/*************************
   SHOW_HELP_INFO
    Prints help info when an invalid command is given or when the
    help flag is given.
**************************/
static void show_help_info(std::string progname) {
  std::cout << "\n"
            << "  C-GLASS usage: " << std::endl;
  std::cout << "    " << progname
            << " params_file --flag1 option1 --flag2 option2 ..." << std::endl;
  std::cout << "\n"
            << "  where flags are one of the following: " << std::endl;
  for (int i = 0; i < n_flags; ++i) {
    std::cout << "    --" << long_options[i].name;
    if (desc[i][1].compare("none") != 0)
      std::cout << " " << desc[i][1];
    std::cout << ", -" << (char)long_options[i].val;
    if (desc[i][1].compare("none") != 0)
      std::cout << " " << desc[i][1];
    std::cout << ", " << desc[i][0] << "\n";
  }
}

static void show_version_info() {
  std::cout << "C-GLASS version no. " << cglass_version << "\n";
}

/*************************
   PARSE_OPTS
    New flag switches should be added below, as well as any information that
    needs to be carried forward in the run_options structure defined above.
    You need to add a case for your flag, as well as the shorthand flag name in
    the string in getopt_long function for the tmp variable
**************************/
static run_options parse_opts(int argc, char *argv[]) {
  if (argc == 1) {
    show_help_info(argv[0]);
    exit(1);
  }

  run_options run_opts;
  int tmp;
  while (1) {
    int option_index = 0;
    tmp = getopt_long(argc, argv, "hvdmaplwbMGcg:r:n:R:", long_options,
                      &option_index);
    if (tmp == -1)
      break;
    switch (tmp) {
    case 0:
      break;
    case 'h':
      show_help_info(argv[0]);
      exit(0);
    case 'v':
      show_version_info();
      exit(0);
    case 'd':
      std::cout << "  Running in debug mode\n";
      run_opts.debug = true;
      break;
    case 'r':
      run_opts.run_name_flag = true;
      run_opts.run_name = optarg;
      break;
    case 'n':
      run_opts.n_run_flag = true;
      run_opts.n_runs = atoi(optarg);
      break;
    case 'g':
      run_opts.graphics_flag = true;
      run_opts.n_graph = atoi(optarg);
      break;
    case 'm':
      run_opts.make_movie = true;
      break;
    case 'G':
      run_opts.auto_graph = true;
      break;
    case 'a':
      run_opts.analysis_flag = true;
      break;
    case 'p':
      run_opts.use_posits = true;
      break;
    case 'l':
      run_opts.load_checkpoint = true;
      break;
    case 'R':
      run_opts.reduce_flag = true;
      run_opts.reduce_factor = atoi(optarg);
      break;
    case 'b':
      run_opts.blank_flag = true;
      break;
    case 'w':
      run_opts.with_reloads = true;
      break;
    case 'M':
      run_opts.single_frame = true;
      run_opts.make_movie = true;
      break;
    case 'c':
      run_opts.convert = true;
      break;
    case '?':
      exit(1);
    default:
      std::cout << "  ERROR: Unrecognized option: " << tmp << "\n";
      exit(1);
    }
  }
  if (optind < argc) {
    run_opts.param_file = argv[optind];
    printf("  Using parameter file: %s\n", argv[optind++]);
  } else {
    printf("ERROR! No parameter file given! Exiting.\n");
    exit(1);
  }
  if (optind < argc) {
    printf("ERROR: Unrecognized context for non-option arguments: ");
    while (optind < argc)
      printf("%s ", argv[optind++]);
    putchar('\n');
    exit(1);
  }
  if (run_opts.make_movie) {
    printf("  Making movies");
    if (run_opts.use_posits) {
      printf(" using posit files");
    }
    if (run_opts.with_reloads) {
      printf(" with available reloads");
    }
    if (run_opts.single_frame) {
      printf(" but generating a single bitmap frame.");
    }
    putchar('\n');
  }
  if (run_opts.analysis_flag) {
    printf("  Running analyses");
    if (run_opts.use_posits) {
      printf(" using posit files");
    }
    if (run_opts.with_reloads) {
      printf(" with available reloads.");
    }
    putchar('\n');
  }
  if (run_opts.reduce_flag) {
    printf("  Reducing output file resolution by a factor of %d\n",
           run_opts.reduce_factor);
  }
  if (run_opts.blank_flag) {
    printf("  Doing a blank run -- generating parameter files without running "
           "simulation.\n");
  }  
  if (run_opts.convert) {
    printf("Converting spec/posit files to text files.\n");
  }
  return run_opts;
}

#endif // _CGLASS_PARSE_FLAGS_H_
