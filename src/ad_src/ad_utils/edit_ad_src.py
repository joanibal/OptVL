"""
    based on autoEditReverse.py and autoEditForward.py in ADflow by G. Kenway
"""

# Import modules
import os
import sys
import re
import argparse


# get input and output directories
parser = argparse.ArgumentParser()
parser.add_argument(
    "mode",
    help="Mode of AD (reverse or forward)",
    type=str,
    choices=["reverse", "forward", "preprocess"],
)
parser.add_argument(
    "--files",
    help="list of source files",
    nargs='+',  # Accept one or more arguments
)
parser.add_argument(
    "--input_dir",
    help="Directory of input source files",
    type=str,
)
parser.add_argument(
    "--output_dir",
    help="Directory of output source files",
    type=str,
)
args = parser.parse_args()

# fmt: off
hand_made_include =  "      INCLUDE 'AVL.INC'\n"\
                     "      INCLUDE 'AVL_ad_seeds.inc'\n"
hand_made_heap_mod = "      use avl_heap_diff_inc\n"
heap_mod =           "      use avl_heap_inc\n"
fake_include =       "      INCLUDE 'AVL_tapenade_fake.inc'\n"
fake_include_diff =       "      INCLUDE 'AVL_tapenade_fake_diff.inc'\n"

# fmt: on


if __name__ == "__main__":
    if args.mode in ["reverse", "forward"]:
        print("Directory of input source files  :", args.input_dir)
        print("Directory of output source files :", args.output_dir)

        
        if args.mode == "reverse":
            # Specify file extension
            file_ext = "_b.f"
            tapenade_include = "INCLUDE 'AVL_b.inc'"
            fake_diff_include = "INCLUDE 'AVL_tapenade_fake_b.inc'"
        elif args.mode == "forward":
            # Specify file extension
            tapenade_include = "INCLUDE 'AVL_d.inc'"
            fake_diff_include = "INCLUDE 'AVL_tapenade_fake_d.inc'"
            file_ext = "_d.f"
            
        # get the maximum number of vortices
        def read_nvmax(path):
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line.upper().startswith("PARAMETER"):
                        # assumes form PARAMETER (NVMAX=5000)
                        start = line.find('(')
                        end = line.find(')', start)
                        if start != -1 and end != -1:
                            inside = line[start+1:end]
                            if inside.upper().startswith("NVMAX="):
                                return int(inside.split("=", 1)[1])
            return None

        # example
        nvmax = read_nvmax("../includes/ADIMEN.INC")

            
        ad_loop_fixes = {
            # only loop over nvor vortices
            'ii1=1,nvmax': 'ii1=1,nvor',
            'ii2=1,nvmax': 'ii2=1,nvor',
            f'ii1=1,{nvmax}': 'ii1=1,nvor',
            f'ii2=1,{nvmax}': 'ii2=1,nvor',
            'ii1=1,ncdim': 'ii1=1,nc',
            'ii2=1,ncdim': 'ii2=1,nc',
            
            # only loop over the surfaces we have
            'ii1=1,nfmax': 'ii1=1,NSURF',
            'ii1=nfmax,1,-1': 'ii1=NSURF,1,-1',
            'ii2=1,nfmax': 'ii2=1,NSURF',
            
            # loop over only the sections we have 
            'ii2=1,nsecmax': 'ii2=1,NSEC(ii1)',
            
            # only loop over the strips we have
            'ii1=1,nsmax': 'ii1=1,NSTRIP',
            'ii2=1,nsmax': 'ii2=1,NSTRIP',
            
            
            
            
        }
        ad_commented_lines = {
            'aoper_b.f': {
            'CALL BUILD_AIC()': 'no need to build the AIC again because we assume analysis is run before',
            'CALL PUSHREAL8ARRAY(wc_gam, 3*nvmax**2)': 'wc_gam is unchanged',
            'CALL POPREAL8ARRAY(wc_gam, 3*nvmax**2)': 'wc_gam unchaged ',
            },
            'asetup_b.f' : {
            '      CALL VVOR(betm, iysym, ysym, izsym, zsym, vrcorec, vrcorew, nvor, ' : 'assume analysis run before',
            '     +          rv1, rv2, lvcomp, chordv, nvor, rc, lvcomp, .false., ' : ' ',
            '     +          wc_gam, nvmax)' : ' ',
            
            }
        }

        
        # modify the AD'd source
        # remove calls to bad includes
        
        for f in os.listdir(args.input_dir):
            if f.endswith(file_ext):
                with open(os.path.join(args.input_dir, f), "r") as fid_org, open(
                    os.path.join(args.output_dir, f), "w"
                ) as fid_mod:
                    print("\nParsing input file", fid_org.name)
                    src_name = os.path.basename(fid_org.name)


                    lines = fid_org.readlines()
                    i = 0
                    while i < len(lines):
                        line = lines[i]
                    
                    
                        # check for bounds fixes
                        for bad_bounds in ad_loop_fixes:
                            if bad_bounds in line:
                                line = line.replace(bad_bounds, ad_loop_fixes[bad_bounds])
                                break
                            
                        # check for commented lines fixes
                        if src_name in ad_commented_lines:
                            comment_dict = ad_commented_lines[src_name]
                            for bad_line in comment_dict:
                                if bad_line in line:
                                    line = f'c     {line[:-1]} {comment_dict[bad_line]}\n'
                                    break
                                
                        if tapenade_include in line:
                            # check to see if the next line is the fake include
                            if fake_diff_include.strip() in lines[i+1]:
                                fid_mod.write(heap_mod)
                                fid_mod.write(hand_made_heap_mod)                                
                                fid_mod.write(hand_made_include)
                                i += 2
                                continue
                            else: 
                                # Insert the hand-written include file
                                fid_mod.write(hand_made_include)

                        elif "REAL*(" in line:
                            # This syntax is not supported by our compiler
                            line = line.replace("REAL*(", "REAL(")
                            fid_mod.write(line)
                        elif "IMPLICIT NONE" in line:
                            # remove the statement
                            pass
                        else:
                            fid_mod.write(line)
                        
                        i += 1

                                
                        # if heap_mod in line:
                        #     # Insert the hand-written diff mod_file
                        #     fid_mod.write(line)
                        #     fid_mod.write(hand_made_heap_mod)
                        # elif fake_diff_include.strip() in line:
                            
                        #     # use the module instead
                        #     fid_mod.write(heap_mod)
                        #     fid_mod.write(hand_made_heap_mod)                        

                print(" Modified file saved", fid_mod.name)

        
        
    if args.mode == "preprocess":
        
        
        # modify the AD'd source
        # remove calls to bad includes
        
        for f in args.files:
            with open(f, "r") as fid_org, open(os.path.join(args.output_dir, os.path.basename(f)), "w") as fid_mod:
                print("\nParsing input file", fid_org.name)

                for line in fid_org:
                    if heap_mod.strip() in line:
                        # use fake include file
                        fid_mod.write(fake_include)
                    elif hand_made_heap_mod.strip() in line:
                        # use fake include file
                        fid_mod.write(fake_include_diff)
                    else:
                        fid_mod.write(line)
                    
    
            print(" Modified file saved", fid_mod.name)

        