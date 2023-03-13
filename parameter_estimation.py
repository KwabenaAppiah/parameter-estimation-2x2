import sys
from linear_nudging_alg import LinearNudgingAlg

def main ():
    if len(sys.argv) == 7:
            str_inputs = sys.argv[1:3]
            ev_type, pp_type = str_inputs
            int_inputs = sys.argv[3:7]
            mu_val, relax_time, bound_val, loop_limit = int(int_inputs[0]), int(int_inputs[1]), int(int_inputs[2]), int(int_inputs[3])
            lin_nudge_alg = LinearNudgingAlg(mu_val, relax_time, ev_type, pp_type, bound_val, loop_limit)
            # bound_val = 5 recommended for testing purposes
            # Note: The algorithm starts to perform poorly w/ a higher bound_val.


    elif len(sys.argv) != 7:
        print("ERROR: You're missing one or more parameters.")

    else:
        print("ERROR: Please check your parameters.")


if __name__ == '__main__':
    main()
