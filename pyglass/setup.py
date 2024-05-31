import os


def make_output_directory(output_directory):
    os.system(f"rm -r {output_directory}")
    os.system(f"mkdir {output_directory}")
    os.system(f"mkdir {output_directory}/ply_files")
    os.system(f"mkdir {output_directory}/temp_images")
    # os.system(f"mkdir {output_directory}/logs")
    # os.system(f"mkdir {output_directory}/config")
    # os.system(f"mkdir {output_directory}/checkpoints")


# def run_cglass(output_directory):
#     cglass_exe_path = f"../output"

# specify the target directory and the base directory
# target_dir = "/path/to/target/directory"
# base_dir = "/path/to/base/directory"

# get the relative path

# %%
os.system("ls")
# make_output_directory(output_directory)
# %%
# all input paths should be relative to C-GLASS/ directory
cglass_dir = "/home/wlough/git/C-GLASS"
os.chdir(cglass_dir)
output_directory = "./output"
params_path = "./meshbrane_tests/mbtest.yaml"
aux_param_paths = [
    "./meshbrane_tests/file1.yaml",
    "./meshbrane_tests/file2.yaml",
    "./meshbrane_tests/cut7_bind.yaml",
]
cglass_exe_path = "./build/src/executable/cglass.exe"
cglass_flags = ""
mkdir_output = True
###########
if mkdir_output:
    make_output_directory(output_directory)
###########
cglass_relative_path = os.path.relpath(cglass_exe_path, output_directory)
params_filename = os.path.basename(params_path)
###########
cp_command = f"cp {params_path} "
for aux_param in aux_param_paths:
    cp_command += f"{aux_param} "
cp_command += output_directory
###########
os.system(cp_command)
os.chdir(output_directory)
run_command = f"{cglass_relative_path} {cglass_flags} {params_filename}"
os.system(run_command)
# %%

import yaml

with open(sim_params_path, "r") as file:
    data = yaml.safe_load(file)

# Accessing a key value
run_name = data["run_name"]
print(run_name)


# %%
import subprocess

# specify the directory you want to change to
dir_path = "./output"

# specify the command you want to run
command = ["ls", "-l"]  # replace with your command

# run the command in the specified directory
result = subprocess.run(command, cwd=dir_path)

# print the return code (0 means success)
print(result.returncode)
