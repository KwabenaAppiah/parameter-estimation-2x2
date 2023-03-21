import sys
from linear_nudging_alg import LinearNudgingAlg  # lorpy
import numpy as np

def main():
    """
    Command line arguments:
      1: String:  Eigenvalue type,
      2: String:  Phase portrait type,
      3: Integer: Mu,
      4: Integer: Relaxation time,
      5: Integer: Bound values,
      6: Integer: Loop limit
    """
    if len(sys.argv) == 7:
        ev_type = sys.argv[1]
        pp_type = sys.argv[2]
        mu_val = int(sys.argv[3])
        relax_time = int(sys.argv[4])
        bound_val = int(sys.argv[5])
        loop_limit = int(sys.argv[6])
        lin_nudge_alg = LinearNudgingAlg(ev_type, pp_type, mu_val, relax_time, bound_val, loop_limit)

    else:
        print("This process requires exactly 7 parameters.")


if __name__ == "__main__":
    main()
