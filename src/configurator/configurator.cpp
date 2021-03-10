#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <yaml-cpp/yaml.h>

int parse_error() {
  fprintf(stderr, "You must provide a valid YAML file to the configurator!\n");
  return 1;
}

class Configurator {
private:
  YAML::Node node_;
  YAML::const_iterator subnode_;

  std::string pfile_;
  std::string pname_;
  std::string pvalue_;
  std::string ptype_;
  std::string subname_;
  std::string whitespace_;
  std::vector<std::string> pfiles_;

  std::ofstream parse_params_;
  std::ofstream parameters_;
  std::ofstream default_config_;
  std::ostringstream store_string_;

  void WriteParameters();
  void WriteParam();
  void WriteParseParams();
  void WriteParseSpeciesParams();
  void WriteParseSpeciesBaseParams();
  void WriteParseSystemParams();
  void WriteDefParam(std::string parent, std::string subparent);
  void WriteDefaultParams();
  void ParseParam(std::string iterator, bool store = false);

public:
  Configurator(std::string default_config_file)
      : parse_params_("include/cglass/parse_params.hpp", std::ios_base::out),
        parameters_("include/cglass/parameters.hpp", std::ios_base::out),
        default_config_("include/cglass/default_params.hpp",
                        std::ios_base::out) {
    node_ = YAML::LoadFile(default_config_file);
  }

  void ConfigureCGlass();
};

struct stat info;

int main(int argc, char *argv[]) {
  if (argc == 1) {
    return parse_error();
  }

  char pathname[] = "include/cglass";
  if (stat(pathname, &info) != 0) {
    printf("Cannot find %s: Are you in the project's root directory?\n", pathname);
    exit(1);
  } else if (info.st_mode & S_IFDIR) {
    printf("Configuring cglass parameters\n");
  } else {
    printf("%s is not a directory.\n", pathname);
    exit(1);
  }
  try {
    Configurator config(argv[1]);
    config.ConfigureCGlass();
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return parse_error();
  }
  return 0;
}

void Configurator::ConfigureCGlass() {
  WriteParameters();
  WriteParseParams();
  WriteDefaultParams();
}

void Configurator::WriteParameters() {
  // Write header
  parameters_ << "#ifndef _CGLASS_PARAMETERS_H_\n"
                 "#define _CGLASS_PARAMETERS_H_\n\n"
                 "#include \"definitions.hpp\"\n\n"
                 "#include <string>\n\n"
                 "template <unsigned char S> struct species_parameters"
                 " {\n";
  // First parse sub parameters (species params, etc)
  for (YAML::const_iterator it = node_.begin(); it != node_.end(); ++it) {
    if (it->first.as<std::string>().compare("species") == 0) {
      if (!it->second.IsMap()) {
        fprintf(stderr,
                "ERROR! Species node is not a yaml map in configurator\n");
      }
      for (subnode_ = it->second.begin(); subnode_ != it->second.end();
           ++subnode_) {
        WriteParam();
      }
    }
  }
  parameters_ << "  virtual ~species_parameters() {}\n"
                 "};\n\n"
                 "typedef species_parameters<species_id::none> species_base"
                 "_parameters;\n\n";

  for (auto it = node_.begin(); it != node_.end(); ++it) {
    std::string name = it->first.as<std::string>();
    if (it->second.IsMap() && name.compare("species") != 0) {
      parameters_ << "template <>\n"
                     "struct species_parameters<species_id::"
                  << name
                  << ">\n"
                     "    : public species_base_parameters {\n";
      for (subnode_ = it->second.begin(); subnode_ != it->second.end();
           ++subnode_) {
        WriteParam();
      }
      parameters_ << "};\n"
                     "typedef species_parameters<species_id::"
                  << name << "> " << name << "_parameters;\n\n";
    }
  }

  // Then parse remaining system parameters
  for (auto it = node_.begin(); it != node_.end(); ++it) {
    parameters_ << "struct system_parameters {\n";
    for (subnode_ = node_.begin(); subnode_ != node_.end(); ++subnode_) {
      if (subnode_->second.IsSequence() && subnode_->second.size() == 2) {
        WriteParam();
      }
    }
    // Write end of file
    parameters_ << "};\n\n#endif // _CGLASS_PARAMETERS_H_";
    parameters_.close();
  }
}

void Configurator::WriteParam() {
  if (!subnode_->second.IsSequence() && subnode_->second.size() != 2) {
    if (subnode_->first.as<std::string>().compare("anchor") == 0) {
      parameters_ << "  struct anchor_parameters {\n";

      YAML::const_iterator subnode_temp = subnode_;

      // Recursively write out anchor parameters
      for (subnode_ = subnode_temp->second.begin(); subnode_ != subnode_temp->second.end(); ++subnode_) {
        parameters_ << "  ";
        WriteParam();
      }
      subnode_ = subnode_temp;
      parameters_ << "  };\n"
                     "  anchor_parameters anchor_temp;\n"
                     "  std::vector<anchor_parameters> anchors = {anchor_temp, anchor_temp};\n";
      return;
    } else {
      fprintf(stderr, "ERROR! Parameter config received an invalid parameter "
                      "that was expected to be sequence of [value, type]\n");
      exit(1);
    }
  }
  pname_ = subnode_->first.as<std::string>();
  pvalue_ = subnode_->second[0].as<std::string>();
  ptype_ = subnode_->second[1].as<std::string>();
  if (ptype_.compare("string") == 0) {
    parameters_ << "  std::string " << pname_ << " = \"" << pvalue_ << "\";\n";
  }
  else if (ptype_.compare("int") == 0 || ptype_.compare("double") == 0 ||
             ptype_.compare("long") == 0 || ptype_.compare("bool") == 0) {
    parameters_ << "  " << ptype_ << " " << pname_ << " = " << pvalue_ << ";\n";
  } else {
    fprintf(stderr, "ERROR! Parameter type %s not recognized!", ptype_.c_str());
    exit(1);
  }
}

void Configurator::WriteDefParam(std::string parent = "", std::string subparent = "") {
  pname_ = subnode_->first.as<std::string>();
  pvalue_ = subnode_->second[0].as<std::string>();
  if (parent.empty()) {
    default_config_ << "  default_config[\"" << pname_ << "\"] = \"" << pvalue_
                    << "\";\n";
  } else {
    if (subparent.empty()) {
      default_config_ << "  default_config[\"" << parent << "\"][\"" << pname_
                      << "\"] = \"" << pvalue_ << "\";\n";
    } else {
      default_config_ << "  default_config[\"" << parent << "\"][\"" << subparent
                      << "\"][\"" << pname_ << "\"] = \"" << pvalue_ << "\";\n";
    }
  }
}

void Configurator::WriteDefaultParams() {
  default_config_ << "  YAML::Node default_config;\n";
  for (YAML::const_iterator it_outer = node_.begin(); it_outer != node_.end(); ++it_outer) {
    if (it_outer->second.IsMap()) {
      std::string parent = it_outer->first.as<std::string>();
      for (YAML::const_iterator it_inner = it_outer->second.begin(); it_inner != it_outer->second.end();
           ++it_inner) {
        if (it_inner->second.IsMap()) {
          std::string subparent = it_inner->first.as<std::string>();
          // Parse subspecies parameters
          for (subnode_ = it_inner->second.begin(); subnode_ != it_inner->second.end(); ++subnode_) {
            WriteDefParam(parent, subparent);
          }
        }
      }
      // Parse species parameters
      for (subnode_ = it_outer->second.begin(); subnode_ != it_outer->second.end(); ++subnode_) {
        if (subnode_->second.IsSequence() && subnode_->second.size() == 2) {
          WriteDefParam(parent, "");
        }
      }
    }
  }
  // Then parse remaining system parameters
  for (subnode_ = node_.begin(); subnode_ != node_.end(); ++subnode_) {
    if (subnode_->second.IsSequence() && subnode_->second.size() == 2) {
      WriteDefParam("", "");
    }
  }
  default_config_.close();
}

void Configurator::WriteParseParams() {
  // Write header
  parse_params_
      << "#ifndef _CGLASS_PARSE_PARAMS_H_\n#define "
         "_CGLASS_PARSE_PARAMS_H_\n\n#include \"yaml-cpp/yaml.h\"\n#include "
         "\"auxiliary.hpp\"\n\n";

  // Write system parameters parser
  WriteParseSystemParams();
  WriteParseSpeciesBaseParams();
  WriteParseSpeciesParams();

  // Write end of file
  parse_params_ << "#endif // _CGLASS_PARSE_PARAMS_H_";
  parse_params_.close();
}

void Configurator::WriteParseSpeciesBaseParams() {
  parse_params_
      << "void parse_species_base_params(species_base_parameters &params,\n"
         "                               YAML::Node &node) {\n"
         "  for (auto it = node.begin(); it != node.end(); ++it) {\n"
         "    if (it->first.as<std::string>().compare(\"species\") == 0) {\n"
         "      if (!it->second.IsMap()) {\n"
         "        Logger::Error(\"Species base params yaml node should always"
         " be a map!\");\n"
         "      }\n"
         "      for (auto jt = it->second.begin(); jt != it->second.end(); "
         "++jt"
         ") {\n"
         "        std::string param_name = jt->first.as<std::string>();\n"
         "        if (false) {\n";
  whitespace_ = "        ";
  for (YAML::const_iterator it = node_.begin(); it != node_.end(); ++it) {
    if (it->first.as<std::string>().compare("species") == 0) {
      if (!it->second.IsMap()) {
        fprintf(stderr,
                "ERROR! Species node is not a yaml map in configurator\n");
        exit(1);
      }
      for (subnode_ = it->second.begin(); subnode_ != it->second.end();
           ++subnode_) {
        if (subnode_->second.IsSequence() && subnode_->second.size() == 2) {
          ParseParam("jt", true);
        }
      }
    }
  }
  parse_params_ << "        } else {\n"
                   "          Logger::Warning(\"Species base parameter %s not"
                   " recognized!\", param_name.c_str());\n"
                   "        }\n"
                   "      }\n"
                   "      return;\n"
                   "    }\n"
                   "  }\n"
                   "}\n\n";
}

void Configurator::WriteParseSpeciesParams() {
  parse_params_ << "species_base_parameters *parse_species_params("
                   "std::string sid,\n"
                   "                                              "
                   "YAML::Node &subnode,\n"
                   "                                              "
                   "YAML::Node &node) {\n  if (false) {\n";
  for (YAML::const_iterator it = node_.begin(); it != node_.end(); ++it) {
    if (it->second.IsMap()) {
      std::string species_name = it->first.as<std::string>();
      if (species_name.compare("species") == 0) {
        continue;
      }
      parse_params_
          << "  } else if (sid.compare(\"" << species_name
          << "\") == 0) {\n    " << species_name
          << "_parameters params;\n"
             "    parse_species_base_params(params, node);\n"
             "    for (auto jt = subnode.begin(); jt != subnode.end();"
             " ++jt) {\n"
             "      std::string param_name = jt->first.as<std::string>();\n"
             "      if (false) {\n";
      // First parse base species parameters
      parse_params_ << store_string_.str();
      // Then write species-specific params
      whitespace_ = "      ";
      for (subnode_ = it->second.begin(); subnode_ != it->second.end();
           ++subnode_) {
        ParseParam("jt");
      }
      parse_params_
          << "      } else {\n"
             "        Logger::Warning(\"Unrecognized %s parameter: '%s'\", "
             "sid.c_str(), param_name.c_str());\n"
             "      }\n"
             "    }\n";
      parse_params_ << "    return new " << species_name
                    << "_parameters(params);\n";
    }
  }
  parse_params_
      << "  } else {\n    Logger::Error(\"Unrecognized SID '%s' in parse_"
         "params!\", sid.c_str());\n  }\n  return nullptr;\n}\n\n";
}

void Configurator::WriteParseSystemParams() {
  parse_params_
      << "system_parameters parse_system_params(YAML::Node "
         "&node) {\n  system_parameters params;\n  for (auto it=node.begin()"
         "; it!= node.end(); ++it) {\n    if (!it->second.IsScalar()) {\n    "
         " "
         "continue;\n    }\n    std::string param_name = "
         "it->first.as<std::string>();\n    if (false) {\n";
  whitespace_ = "    ";
  for (subnode_ = node_.begin(); subnode_ != node_.end(); ++subnode_) {
    if (subnode_->second.IsSequence() && subnode_->second.size() == 2) {
      ParseParam("it");
    }
  }
  parse_params_ << "    } else {\n      Logger::Warning(\"Unrecognized "
                   "parameter '%s'\", param_name.c_str());\n    }\n  }\n"
                   "  return params;\n}\n\n";
}

void Configurator::ParseParam(std::string iterator, bool store) {
  pname_ = subnode_->first.as<std::string>();
  if (pname_.compare("anchor") == 0) {
    parse_params_ << "      } else if (param_name.compare(\"" << pname_ << "\")==0) {\n"
                     "        int index = 0;\n"
                     "        for (auto seq=jt->second.begin(); seq!= jt->second.end(); ++seq) {\n"
                     "          if (index > 1) {\n"
                     "            Logger::Error(\"Only two anchors allowed per crosslink.\");\n"
                     "          }\n"
                     "          for (auto kt = seq->second.begin(); kt != seq->second.end(); ++kt) {\n"
                     "            if (!kt->second.IsScalar())  {\n"
                     "              continue;\n"
                     "            }\n"
                     "            std::string sub_param_name = kt->first.as<std::string>();\n"
                     "            if (false) {\n";
    for (YAML::const_iterator subit = subnode_->second.begin(); subit != subnode_->second.end();
           ++subit) {
      if (subit->second.IsSequence() && subit->second.size() == 2) {
        pname_ = subit->first.as<std::string>();
        ptype_ = subit->second[1].as<std::string>();
        if (ptype_.compare("string") == 0) {
          ptype_ = "std::string";
        }
        parse_params_ <<
                         "            } else if (sub_param_name.compare(\"" << pname_ << "\")==0) {\n"
                         "              params.anchors[index]." << pname_ << " = kt->second.as<"
                         << ptype_ << ">();\n";
      }
    }
    parse_params_ << "            } else {\n              "
                     "Logger::Warning(\"Unrecognized parameter '%s'\", sub_param_name.c_str());\n"
                     "            }\n          }\n          index++;\n        }\n";
    return;
  }
  ptype_ = subnode_->second[1].as<std::string>();
  if (ptype_.compare("string") == 0) {
    ptype_ = "std::string";
  }
  parse_params_ << whitespace_ << "} else if (param_name.compare(\"" << pname_
                << "\")==0) {\n"
                << whitespace_ << "params." << pname_ << " = " << iterator
                << "->second.as<" << ptype_ << ">();\n";
  if (store) {
    store_string_ << "      } else if (param_name.compare(\"" << pname_
                  << "\")==0) {\n"
                  << "      params." << pname_ << " = " << iterator
                  << "->second.as<" << ptype_ << ">();\n";
  }
}
