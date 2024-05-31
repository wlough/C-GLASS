# %%
import os
import shutil
import subprocess
import numpy as np
import yaml


# %%


class SimWrapper:
    def __init__(
        self,
        cglass_params_path,
        cglass_insert_paths,
        name=None,
        C_GLASS_ROOT_DIR=None,
        C_GLASS_OUTPUT_DIR=None,
        regenerate_output_dir=False,
    ):
        self.cglass_params_path = cglass_params_path
        self.cglass_insert_paths = cglass_insert_paths

        self.cglass_spec_keys_all = ["spindle", "filament", "crosslink", "brane"]
        if C_GLASS_ROOT_DIR is None:
            self.C_GLASS_ROOT_DIR = "/home/wlough/git/C-GLASS"
        else:
            self.C_GLASS_ROOT_DIR = C_GLASS_ROOT_DIR
        if C_GLASS_OUTPUT_DIR is None:
            self.C_GLASS_OUTPUT_DIR = os.path.join(self.C_GLASS_ROOT_DIR, "output")
        else:
            self.C_GLASS_OUTPUT_DIR = C_GLASS_OUTPUT_DIR

        # self.cglass_exe_path = os.path.join(C_GLASS_ROOT_DIR, "build", "src", "executable", "cglass.exe")

        with open(cglass_params_path, "r") as f:
            self.cglass_params = yaml.load(f, Loader=yaml.FullLoader)

        if name is None:
            self.name = self.cglass_params["run_name"]
        else:
            self.name = name

        ######################################
        self.check_output_dir(regenerate=regenerate_output_dir)
        self.copy_yamls()

        self.cglass_insert_files = {}
        for insert_path in cglass_insert_paths:
            with open(insert_path, "r") as f:
                insert_file_name = os.path.basename(insert_path)
                self.cglass_insert_files[insert_file_name] = yaml.load(f, Loader=yaml.FullLoader)

        self.spec_keys = [key for key in self.cglass_spec_keys_all if key in self.cglass_params.keys()]

    def check_output_dir(self, regenerate=False):
        """checks/generates the output directory and subdirectories for the simulation."""
        C_GLASS_OUTPUT_DIR = self.C_GLASS_OUTPUT_DIR
        output_dir = self.get_output_dir()
        subdirs = ["ply_files", "temp_images", "checkpoints"]
        if regenerate and os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(C_GLASS_OUTPUT_DIR, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        for subdir in subdirs:
            _subdir = os.path.join(output_dir, subdir)
            os.makedirs(_subdir, exist_ok=True)

    def copy_yamls(self):
        file_paths = [self.cglass_params_path, *self.cglass_insert_paths]
        target_dir = self.get_output_dir()
        for f in file_paths:
            if os.path.dirname(f) != target_dir:
                # target_file_path = os.path.join(target_dir, os.path.basename(f))
                # if not os.path.exists(target_file_path):
                shutil.copy(f, target_dir)

    def get_output_dir(self):
        return os.path.join(self.C_GLASS_OUTPUT_DIR, self.name)

    def get_cglass_exe_path(self):
        return os.path.join(self.C_GLASS_ROOT_DIR, "build", "src", "executable", "cglass.exe")

    def run_cglass_exe(self, cglass_flags=[]):

        # cwd = os.getcwd()
        cglass_exe = self.get_cglass_exe_path()
        output_dir = self.get_output_dir()
        params = os.path.basename(self.cglass_params_path)
        # os.chdir(output_dir)
        # run_command = f"{cglass_exe_path} {cglass_flags} {params_name}"
        # subprocess.run([cglass_exe_path, cglass_flags, params_name])
        # Start the process
        # dir_path = "path/to/your/directory"
        run_command = [cglass_exe, params] + cglass_flags

        # Start the process
        process = subprocess.Popen(
            run_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=output_dir,
        )
        # process = subprocess.Popen(exe_path, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        # Print the output in real-time
        for line in iter(process.stdout.readline, b""):
            print(line.decode(), end="")


# %%

cglass_params_path = "/home/wlough/git/C-GLASS/examples/spindle_params.yaml"  # "./pyglass/params.yaml"
cglass_insert_paths = ["/home/wlough/git/C-GLASS/examples/assay_3d_insert.yaml"]
# cglass_params_path = "/home/wlough/git/C-GLASS/examples/cut7_crosslinking.yaml"  # "./pyglass/params.yaml"
# cglass_insert_paths = [
#     "/home/wlough/git/C-GLASS/examples/cut7_bind.yaml",
#     "/home/wlough/git/C-GLASS/examples/file1.yaml",
#     "/home/wlough/git/C-GLASS/examples/file2.yaml",
# ]

s = SimWrapper(cglass_params_path, cglass_insert_paths, regenerate_output_dir=True)
s.run_cglass_exe()

# %%
