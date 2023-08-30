from linear_nudging_alg import LinearNudgingAlg
import numpy as np
import shutil
import sys
import os


def main():
    delete_output_dir("../output/")

    if len(sys.argv) == 5: #actually 6, but 7 in total
        mu_val = int(round(float(sys.argv[1])))
        relax_time = float(sys.argv[2])
        case_type = sys.argv[3]
        file_import = sys.argv[4]

        substrings = file_import.split('/')
        filename = substrings[2]
        str_list  = filename.split("_")
        ev_type =  str_list[0]
        pp_type = str_list[1]



        if ev_type == "rde":
            if pp_type == "saddle" or pp_type == 'sink' or pp_type == 'source':
                lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, case_type, file_import)
            else:
                print("Please select a phase portrait type of 'saddle', 'sink' or 'source'.")

        elif ev_type == "re":
            if pp_type == 'sink' or pp_type == 'source':
                lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, case_type, file_import)
            else:
                print("Please select a phase portrait type of 'sink' or 'source'.")

        elif ev_type == "ce":

            if str_list[1] == "sp" and str_list[2] == "sink" or str_list[2] == "source" :
                pp_type  = str_list[1] + "_"  + str_list[2]
                
                lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, case_type, file_import)

            elif pp_type == 'center':
                    lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, case_type, file_import)
            else:
                print("Please select a phase portrait type of 'sp_sink', 'sp_source' or 'center'.")

        else:
            print("Please specifiy an eigenvalue type of 'rde', 're' or 'ce'.")


    elif len(sys.argv) == 8: # actually 6, but 7 in total
        ev_type = sys.argv[1]
        pp_type = sys.argv[2]
        mu_val = int(round(float(sys.argv[3])))
        relax_time = float(sys.argv[4])
        bound_val = int(sys.argv[5])
        loop_limit = int(sys.argv[6])
        case_type = str(sys.argv[7])
        lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, case_type, bound_val, loop_limit)

    else:
        print("This application requires either 5 or 7 parameters.")


def delete_output_dir(path):
    # Check if the directory exists
    if os.path.exists(path):
        # Remove the directory
        shutil.rmtree(path)
    else:
        print("")
        print("WARNING: The directory", "'" + path + "'", "does not exist.", "\n")

if __name__ == "__main__":
    main()
