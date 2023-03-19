import sys
from linear_nudging_alg import LinearNudgingAlg  # lorpy
import numpy as np

def main():
    """
    command line arguments:
      1: string:  eigenvalue type,
      2: string:  phase portrait type,
      3: integer: mu,
      4: integer: relaxation time,
      5: integer: bound_val, # <= 5 for best performance
      6: integer: loop_limit
    """
    if len(sys.argv) == 7:
        ev_type = sys.argv[1]
        pp_type = sys.argv[2]
        mu_val = int(sys.argv[3])
        relax_time = int(sys.argv[4])
        bound_val = float(sys.argv[5])
        loop_limit = int(sys.argv[6])
        lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, bound_val, loop_limit)

    else:
        print("This process requires exactly 7 parameters.")


if __name__ == "__main__":
    main()
