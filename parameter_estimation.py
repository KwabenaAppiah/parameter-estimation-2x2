import sys
from linear_nudging_alg import LinearNudgingAlg  # lorpy
import numpy as np


def main():

    # if len(sys.argv) == 6: #actually 5, but 6 in total
    #     import_file_path = sys.argv[5]
    #     lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, import_file_path)

    if len(sys.argv) == 4: #actually 5, but 6 in total
        mu_val = int(round(float(sys.argv[1])))
        relax_time = int(sys.argv[2])
        file_import = sys.argv[3]

        substrings = file_import.split('/')
        filename = substrings[2]
        str_list  = filename.split("_")
        ev_type =  str_list[0]
        pp_type = str_list[1]
    

        if ev_type == "rde":
            if pp_type == "saddle" or pp_type == 'sink' or pp_type == 'source':
                lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, file_import)
            else:
                print("Please select a phase portrait type of 'saddle', 'sink' or 'source'.")

        elif ev_type == "re":
            if pp_type == 'sink' or pp_type == 'source':
                lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, file_import)
            else:
                print("Please select a phase portrait type of 'sink' or 'source'.")

        elif ev_type == "ce":

            if str_list[1] == "sp" and str_list[2] == "sink" or str_list[2] == "source" :
                pp_type  = str_list[1] + "_"  + str_list[2]
                lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, file_import)

            elif pp_type == 'center':
                    lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, file_import)
            else:
                print("Please select a phase portrait type of 'sp_sink', 'sp_source' or 'center'.")

        else:
            print("Please specifiy an eigenvalue type of 'rde', 're' or 'ce'.")



    elif len(sys.argv) == 7: # actually 6, but 7 in total
        ev_type = sys.argv[1]
        pp_type = sys.argv[2]
        mu_val = int(round(float(sys.argv[3])))
        relax_time = int(sys.argv[4])
        bound_val = int(sys.argv[5])
        loop_limit = int(sys.argv[6])
        lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, bound_val, loop_limit)

    else:
        # print("This application requires a minimum of 6 parameters and a max of 7.")
        print("This application requires either 4 or 7 parameters.")




if __name__ == "__main__":
    main()
