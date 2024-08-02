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

    areas_to_cd_into = {}
    for name, area in areas_to_merge.items():
        areas_to_cd_into[name] = str(os.path.abspath(area)) + "/SubProcesses"
    
    coupling_files_to_access = {}
    for name in areas_to_merge.keys():
        coupling_files_to_access[name] = PWD + "/archive/" + name + "_coupl.inc" 
    
    
    coupling_files_to_replace = {}
    for name, area in areas_to_merge.items():
        coupling_files_to_replace[name] = str(os.path.abspath(area)) + "/Source/MODEL/coupl.inc" 
    
    
    input_files_to_access = {}
    for name in areas_to_merge.keys():
        input_files_to_access[name] = PWD + "/archive/" + name + "_input.inc" 
    
    
    input_files_to_replace = {}
    for name, area in areas_to_merge.items():
        input_files_to_replace[name] = str(os.path.abspath(area)) + "/Source/MODEL/input.inc" 
    

    makefile_options_to_replace = {}
    for name, area in areas_to_merge.items():
        makefile_options_to_replace[name] = str(os.path.abspath(area)) + "/Source/make_opts" 
    

    ##################### COMBINING COMMON BLOCKS #################3

    for name in coupling_files_to_replace.keys():
        coupling_file = coupling_files_to_replace[name]
        coupling_file_name = coupling_file.split("/")[-1]
        if not os.path.islink(coupling_file):
            subprocess.run(["cp", coupling_file, "archive/" + name + "_" + coupling_file_name], check=True)
        
        input_file = input_files_to_replace[name]
        input_file_name = input_file.split("/")[-1]
        if not os.path.islink(input_file):
            subprocess.run(["cp", input_file, "archive/" + name + "_" + input_file_name], check=True)


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
            named_couplings = {}
            for match in all_matches:
                named_couplings[match.group("COMMON_NAME").strip()] = "".join(match.group("terms").split()).replace('$', '')

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
    with open(output_area + "/coupl.inc", "w+") as f:
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
            named_inputs = {}
            for match in all_matches:
                named_inputs[match.group("COMMON_NAME").strip()] = "".join(match.group("terms").split()).replace('$', '')

            for input_name, input_set in named_inputs.items():
                input_set = set(input_set.strip().split(","))
                if input_name in inputs.keys():
                    inputs[input_name] |= input_set
                else:
                    inputs[input_name] = input_set

    for name, input_new in inputs.items():
        inputs[name] = sorted(input_new)


    print("Creating combined input.inc...")
    with open(output_area + "/input.inc", "w+") as f:
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

            f.write("\n\n\n" + " "*6 + "COMMON/" + name + "/ ")

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
        subprocess.call(["rm", makefile_options])
        os.symlink("../../libSMEFTSIM/make_opts", makefile_options)

        os.chdir("MODEL")
        coupling_file = coupling_files_to_replace[name]
        subprocess.check_call(["rm", coupling_file])
        os.symlink("../../../libSMEFTSIM/coupl.inc", coupling_file)

        input_file = input_files_to_replace[name]
        subprocess.check_call(["rm", input_file])
        os.symlink("../../../libSMEFTSIM/input.inc", input_file)

        subprocess.check_call("make clean", shell=True)
        subprocess.check_call("make", shell=True)
        os.chdir("../DHELAS")
        subprocess.check_call("make clean", shell=True)
        subprocess.check_call("make", shell=True)


    for name, area in areas_to_cd_into.items():
        print("Compiling", name)

        os.chdir(area)
        subprocess.call("rm *.o */*.o", shell=True)
        
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
        
        subprocess.check_call(["make", "cpp"])
        subprocess.call("rm *.a", shell=True)
        files_to_compile = []
        for compilation_glob in (
            "../Source/DHELAS/*.o", 
            "../Source/MODEL/*.o", 
            "all_matrix.o", 
            "*/matrix.o"):
            files_to_compile += glob.glob(compilation_glob)
        subprocess.check_call(["ar", "cru", "lib" + name + ".a"] + files_to_compile)
        subprocess.check_call(["ranlib", "lib" + name + ".a"])
        p = subprocess.Popen(["nm", "lib" + name + ".a"], stdout=subprocess.PIPE, encoding='utf-8')
        named_values, _ = p.communicate()
        
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
                things_to_rename[columns[2]] = name + "_" + columns[2]
        rename_command = " --redefine-sym "
        temp_command = []
        for old, new in things_to_rename.items():
            temp_command.append(old + "=" + new)
        command = rename_command + rename_command.join(temp_command)
        command = "objcopy lib" + name + ".a lib" + name + "_renamed.a" + command
        # print()
        # print(command)
        # print()
        subprocess.check_call(command, shell=True)
    
    os.chdir(PWD)
    subprocess.call(
        "g++ -Wl,--whole-archive */SubProcesses/*{SIG,BKG,BSI}_renamed.a -Wl,--no-whole-archive -fPIC -lgfortran -shared -o libSMEFTSIM/libSMEFTsim.so > debug_compileAll 2>&1",
        shell=True
    )
    
        
    #g++ -Wl,--whole-archive */SubProcesses/*{SIG,BKG,BSI}_renamed.a -Wl,--no-whole-archive -fPIC -lgfortran -shared -o libSMEFTsim.so > debug_compileAll 2>&1
    # 1) compile each area into a single static library .a file
    # 2) start renaming symbols using objcopy to reflect process numbering scheme
