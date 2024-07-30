import os
import subprocess
import shutil
import sys
import glob
import re
from constants import *

coupling_block_header = """ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION G
      COMMON/STRONG/ G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R
      COMMON/RSCALE/ MU_R

      DOUBLE PRECISION NF
      PARAMETER(NF=0)

      DOUBLE PRECISION MDL_MZ,MDL_MW,MDL_MT,MDL_MB,MDL_MH,MDL_MW1
     $ ,MDL_MZ1,MDL_ME,MDL_MMU,MDL_MTA,MDL_MH1,MDL_MS,MDL_MD,MDL_MU
     $ ,MDL_MC,MDL_MT1

      COMMON/MAD_MASSES/ MDL_MZ,MDL_MW,MDL_MT,MDL_MB,MDL_MH,MDL_MW1
     $ ,MDL_MZ1,MDL_ME,MDL_MMU,MDL_MTA,MDL_MH1,MDL_MS,MDL_MD,MDL_MU
     $ ,MDL_MC,MDL_MT1


      DOUBLE PRECISION MDL_WZ1,MDL_WZ,MDL_WH1,MDL_WT,MDL_WH,MDL_WT1
     $ ,MDL_WW,MDL_WW1

      COMMON/WIDTHS/ MDL_WZ1,MDL_WZ,MDL_WH1,MDL_WT,MDL_WH,MDL_WT1
     $ ,MDL_WW,MDL_WW1



""" #every process shares this part of the coupling header

input_block_header = """ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

""" #every process shares this part of the input header

if __name__ == "__main__":
    PWD = os.path.abspath(os.curdir)
    #This dictates how the processes will be numbered!
    # 0, 1, 2, 3, etc.

    output_area = os.path.abspath("libSMEFTSIM")

    areas_to_cd_into = {
        name : str(os.path.abspath(area)) + "/SubProcesses" for 
        name, area in areas_to_merge.items()
    }
    
    coupling_files_to_access = {
        name : f"{PWD}/archive/{name}_coupl.inc" for 
        name in areas_to_merge.keys()
    }
    coupling_files_to_replace = {
        name : str(os.path.abspath(area)) + "/Source/MODEL/coupl.inc" for 
        name, area in areas_to_merge.items()
    }
    
    input_files_to_access = {
        name : f"{PWD}/archive/{name}_input.inc" for 
        name in areas_to_merge.keys()
    }
    input_files_to_replace = {
        name : str(os.path.abspath(area)) + "/Source/MODEL/input.inc" for 
        name, area in areas_to_merge.items()
    }

    makefile_options_to_replace = {
        name : str(os.path.abspath(area)) + "/Source/make_opts" for
        name, area in areas_to_merge.items()
    }

    ##################### COMBINING COMMON BLOCKS #################3

    for name in coupling_files_to_replace.keys():
        coupling_file = coupling_files_to_replace[name]
        coupling_file_name = coupling_file.split("/")[-1]
        if not os.path.islink(coupling_file):
            subprocess.run(["cp", coupling_file, f"archive/{name}_{coupling_file_name}"], check=True)
        
        input_file = input_files_to_replace[name]
        input_file_name = input_file.split("/")[-1]
        if not os.path.islink(input_file):
            subprocess.run(["cp", input_file, f"archive/{name}_{input_file_name}"], check=True)


    common_block_regex = re.compile(
        r"COMMON\/(?P<COMMON_NAME>\S+)\/(?P<terms>(?:.+\s+)(?:\$\s+.+\s*)*)"
    )


    # start with coupl.inc for each area
    coupling_types = {
        "COUPLINGS":"DOUBLE COMPLEX"
    }
    couplings = {}
    
    for name, coupling_file in coupling_files_to_access.items():
        print("Processing", name)
        
        with open(coupling_file) as f:
            all_matches = re.finditer(common_block_regex, f.read())
            named_couplings = {match.group("COMMON_NAME").strip() :
                "".join(match.group("terms").split()).replace('$', '')
                for match in all_matches
            }
            coupling_name = "COUPLINGS"
            coupling_set = named_couplings[coupling_name]
            coupling_set = set(coupling_set.strip().split(","))
            if coupling_name in couplings.keys():
                couplings[coupling_name] |= coupling_set
            else:
                # print(f"ADDING COUPLINGS TO {coupling_name}", coupling_set)
                couplings[coupling_name] = coupling_set

    couplings["COUPLINGS"] = sorted(
        couplings['COUPLINGS'], key=lambda x: int(x[x.find("_") + 1:]) #sort off number
        )

    print("Creating combined coupl.inc...")
    with open(f"{output_area}/coupl.inc", "w+") as f:
        f.write(coupling_block_header)
        f.write(" "*6 + "DOUBLE COMPLEX ")

        firstLine = True
        i = 0
        for n, coupling in enumerate(couplings["COUPLINGS"]):
            if firstLine and (i % COMMON_BLOCK_FIRST_LINE_LIMIT_TERMS == 0) and i != 0:
                f.write(coupling + ',\n' + " "*5 + '$ ')
                firstLine = False
                i = 0
            elif not firstLine and (i % COMMON_BLOCK_LIMIT_TERMS == 0) and (i != 0) and (n != len(couplings["COUPLINGS"]) - 1):
                f.write(coupling + ',\n' + " "*5 + '$ ')
                i = 0
            else:
                if n == len(couplings["COUPLINGS"]) - 1:
                    f.write(coupling)
                else:
                    f.write(coupling + ',')
            i += 1

        f.write("\n\n\n" + " "*6 + "COMMON/COUPLINGS/ ")

        firstLine = True
        i = 0
        for n, coupling in enumerate(couplings["COUPLINGS"]):
            if firstLine and (i % COMMON_BLOCK_FIRST_LINE_LIMIT_TERMS == 0) and i != 0:
                f.write(coupling + ',\n' + " "*5 + '$ ')
                firstLine = False
                i = 0
            elif not firstLine and (i % COMMON_BLOCK_LIMIT_TERMS == 0) and (i != 0) and (n != len(couplings["COUPLINGS"]) - 1):
                f.write(coupling + ',\n' + " "*5 + '$ ')
                i = 0
            else:
                if n == len(couplings["COUPLINGS"]) - 1:
                    f.write(coupling)
                else:
                    f.write(coupling + ',')
            i += 1
    
    #now let's do input.inc for each area
    input_types = {
        "PARAMS_R":"DOUBLE PRECISION",
        "PARAMS_C":"DOUBLE COMPLEX"
    }
    inputs = {}
    
    for name, input_file in input_files_to_access.items():
        print("Processing", name)
        
        with open(input_file) as f:
            all_matches = re.finditer(common_block_regex, f.read())
            named_inputs = {match.group("COMMON_NAME").strip() :
                "".join(match.group("terms").split()).replace('$', '')
                for match in all_matches
            }
            for input_name, input_set in named_inputs.items():
                input_set = set(input_set.strip().split(","))
                if input_name in inputs.keys():
                    inputs[input_name] |= input_set
                else:
                    inputs[input_name] = input_set

    for name, input_new in inputs.items():
        inputs[name] = sorted(input_new)


    print("Creating combined input.inc...")
    with open(f"{output_area}/input.inc", "w+") as f:
        f.write(input_block_header)
        for name, input_type in input_types.items():
            f.write("\n\n\n" + " "*6 + input_type + " ")

            firstLine = True
            i = 0
            for n, input in enumerate(inputs[name]):
                if firstLine and (i % COMMON_BLOCK_FIRST_LINE_LIMIT_TERMS_input == 0) and i != 0:
                    f.write(input + ',\n' + " "*5 + '$ ')
                    firstLine = False
                    i = 0
                elif not firstLine and (i % COMMON_BLOCK_LIMIT_TERMS_input == 0) and (i != 0) and (n != len(inputs[name]) - 1):
                    f.write(input + ',\n' + " "*5 + '$ ')
                    i = 0
                else:
                    if n == len(inputs[name]) - 1:
                        f.write(input)
                    else:
                        f.write(input + ',')
                i += 1

            f.write("\n\n\n" + " "*6 + f"COMMON/{name}/ ")

            firstLine = True
            i = 0
            for n, input in enumerate(inputs[name]):
                if firstLine and (i % COMMON_BLOCK_FIRST_LINE_LIMIT_TERMS_input == 0) and i != 0:
                    f.write(input + ',\n' + " "*5 + '$ ')
                    firstLine = False
                    i = 0
                elif not firstLine and (i % COMMON_BLOCK_LIMIT_TERMS_input == 0) and (i != 0) and (n != len(inputs[name]) - 1):
                    f.write(input + ',\n' + " "*5 + '$ ')
                    i = 0
                else:
                    if n == len(inputs[name]) - 1:
                        f.write(input)
                    else:
                        f.write(input + ',')
                i += 1

    # make new symbolic links and add the files as references
    print("Creating symbolic links to new files...")
    for name in coupling_files_to_replace.keys():
        makefile_options = makefile_options_to_replace[name]
        os.chdir(makefile_options[:makefile_options.rfind('/')])
        subprocess.run(["rm", makefile_options])
        os.symlink(f"../../libSMEFTSIM/make_opts", makefile_options)

        os.chdir("MODEL")
        coupling_file = coupling_files_to_replace[name]
        subprocess.run(f"rm {coupling_file}", shell=True, check=True)
        os.symlink(f"../../../libSMEFTSIM/coupl.inc", coupling_file)

        input_file = input_files_to_replace[name]
        subprocess.run(f"rm {input_file}", shell=True, check=True)
        os.symlink(f"../../../libSMEFTSIM/input.inc", input_file)

        subprocess.run("make clean", shell=True, check=True)
        subprocess.run("make", shell=True, check=True)
        os.chdir("../DHELAS")
        subprocess.run("make clean", shell=True, check=True)
        subprocess.run("make", shell=True, check=True)


    for name, area in areas_to_cd_into.items():
        print("Compiling", name)

        os.chdir(area)
        subprocess.run("rm *.o */*.o", shell=True)
        
        make_cpp_flag = True
        full_text = ""
        with open("makefile") as makefile:
            for line in makefile:
                if "cpp: $(LIBDIR)/$(LIBRARY)" in line:
                    make_cpp_flag = False #Do not need to make a cpp flag
                full_text += line

        if make_cpp_flag:
            with open("makefile", "a") as makefile:
                makefile.write("\ncpp: $(LIBDIR)/$(LIBRARY)")
        
        subprocess.run(["make", "cpp"], check=True)
        subprocess.run("rm *.a", shell=True)
        files_to_compile = []
        for compilation_glob in (
            "../Source/DHELAS/*.o", 
            "../Source/MODEL/*.o", 
            "all_matrix.o", 
            "*/matrix.o"):
            files_to_compile += glob.glob(compilation_glob)
        subprocess.run(["ar", "cru", f"lib{name}.a"] + files_to_compile, check=True)
        subprocess.run(["ranlib", f"lib{name}.a"], check=True)
        named_values = subprocess.check_output(["nm", f"lib{name}.a"]).decode(sys.stdout.encoding).strip()
        
        all_entries = named_values.split('\n')
        things_to_rename = {}
        for entry in all_entries:
            columns = entry.strip().split()
            if len(columns) == 3 and columns[1] in ("T", "C", "D", "R", "A", "B"):
                if columns[2] in ( # these things are common amongst packages
                    "couplings_", "mad_masses_", 
                    "widths_", "masses_", "strong_", "weak_", "rscale_",
                    "params_r_", "params_c_"
                    ):
                    continue
                things_to_rename[columns[2]] = f"{name}_{columns[2]}"
        rename_command = " --redefine-sym "
        command = rename_command + rename_command.join(
            [f"{old}={new}" for old, new in things_to_rename.items()]
            )
        command = f"objcopy lib{name}.a lib{name}_renamed.a" + command
        # print()
        # print(command)
        # print()
        subprocess.run(command, shell=True, check=True)
    
    os.chdir(PWD)
    subprocess.run(
        "g++ -Wl,--whole-archive */SubProcesses/*{SIG,BKG,BSI}_renamed.a -Wl,--no-whole-archive -fPIC -lgfortran -shared -o libSMEFTSIM/libSMEFTsim.so > debug_compileAll 2>&1",
        shell=True,
        check=True
    )
    
        
    #g++ -Wl,--whole-archive */SubProcesses/*{SIG,BKG,BSI}_renamed.a -Wl,--no-whole-archive -fPIC -lgfortran -shared -o libSMEFTsim.so > debug_compileAll 2>&1
    # 1) compile each area into a single static library .a file
    # 2) start renaming symbols using objcopy to reflect process numbering scheme
