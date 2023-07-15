# parameter-estimation

BACKGROUND

Based on the findings of Elizabeth Carlson, Joshua Hudson and Adam Larios (https://epubs.siam.org/doi/10.1137/19M1248583),
this project utilizes a nudging algorithm for the purpose of parameter recovery. As of now, this script only works for two-dimensional linear systems.


--------------------------------------------------------------------------------------------------


HOW DO I RUN THIS SCRIPT?

There are two command line prompt format options:

Option 1: 
python parameter_estimation.py EV_TYPE PP_TYPE MU_VAL RLX_TIME BND_VAL LOOP_LIM
python parameter_estimation.py ce sp_sink 10000 10 12 200 (Sample prompt)

Option 2: 
python parameter_estimation.py MU_VAL RLX_TIME FILE_IMPORT
python parameter_estimation.py 100 10 ../import/ce_sp_sink_all_matrices.txt
* NOTE: There must be a file within the 'import' directory in order to run this command.



--------------------------------------------------------------------------------------------------


WHAT DOES THIS SCRIPT OUTPUT?

- TBD 


--------------------------------------------------------------------------------------------------
Updated 2023.07.15 - 01:00
