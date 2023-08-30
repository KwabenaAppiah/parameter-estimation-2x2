# parameter-estimation

 BACKGROUND

Based on the findings of Elizabeth Carlson, Joshua Hudson and Adam Larios  (https://epubs.siam.org/doi/10.1137/19M1248583),
this project utilizes a nudging algorithm for the purpose of parameter recovery in a 2 x 2 linear system.

--------------------------------------------------------------------------------------------------

HYPOTHESIS

The nudging algorithm’s performance in recovering parameters in a given matrix are dependent on the dynamics of said matrix's solution.

--------------------------------------------------------------------------------------------------

TESTING

In order to test this hypothesis, the nudging algorithm is applied to an unbiased sample of random matrices within the trace-determinant plane. In each case, this sample of matrices is pulled based on eight possible dynamical behaviors. 

--------------------------------------------------------------------------------------------------

HOW DO I RUN THIS SCRIPT?

There are two command line prompt options:

Option #1 - Generate a random sample of matrices:  

python parameter_estimation.py EV_TYPE PP_TYPE MU_VAL RLX_TIME BND_VAL LOOP_LIM CASE_TYPE 

(E.g. python parameter_estimation.py ce sp_sink 1000 .5 12 5 main_diagonal)

Option #2 - Pull from an existing list of matrices: 

python parameter_estimation.py MU_VAL RLX_TIME FILE_IMPORT 

(E.g. python parameter_estimation.py 1000 .5 main_diagonal ../import/ce_sp_sink_all_matrices.txt)

* NOTE: 
 - There must be a file within the "import" directory in order to run option #2
- Each time option #1 or #2 is run, a text file that can be used for either option, is generated in the directory "text_files".


EV_TYPE (eigenvalue type) and PP_TYPE (phase portrait type):

1. rde = Real distinct eigenvalues
- saddle
- source
- sink

2. re = Repeated eigenvalues
- sink
- source

3. ce =  Complex eigenvalues
- sp_sink = Spiral sink, 
- sp_source = Spiral source
- center

MU_VAL (µ)
- The nudging parameter

RLX_TIME: 
- The relaxation time or relaxation parameter is a static value that is used to determine when to call the update formula within the nudging algorithm.
- When the difference between the current time, t, and the last update time, t_n, elapses the relax time, (i.e. t - t_n => RLX_TIME), the update portion of the nudging algorithm is called.

BND_VAL: 
- The bound value represents the boundaries in which the randomly generated matrix values are pulled from.
- Thus, in the sample line above, we utilize a boundary between -12 and 12, inclusive. 

LOOP_LIM: 
- The total number of randomly generated matrices in which the nudging algorithm will be applied to.

FILE_IMPORT
- An imported text file that contains a list of matrices and their average relative parameter error.

CASE_TYPE
- Given a 2 x 2 matrix, the case type refers to the specific row of parameters we'd like to recover.
- Currently, this application can recover the main diagonal, anti-diagonal, left column and right column of a 2 x 2 matrix.
- Moreover, the following terms should be used in the command line prompts, respectively: main_diagonal, anti-diagonal, left_column and right_column.


--------------------------------------------------------------------------------------------------


WHAT DOES THIS SCRIPT OUTPUT?

Individual Graphs:
- The followings directories provide graphs for each individual matrix. 
- In addition, subfolders within this directory may appended with the words "good" or "bad".
- In this case, “good” indicates a certain matrix is above a certain threshold value, while “bad” signifies said matrix falls below this value. 

1. avg_param_err_graphs 
- Graphs in this directory display the combined average parameter error for each missing parameter.

2. param_err_graphs
- This directory contains graphs in which both the relative and absolute parameter error for each missing parameter is displayed.

3. sol_signal_err_graphs
- Graphs in this directory display the signal error for both the assimilated and true matrix parameters are graphed in one plot.

4. sol_signal_err_split_graphs
- This directory contains graphs in which the signal error for the assimilated and true solution are graphed in two separate plots.

5. sol_x_graphs
- Graphs in this directory display the x value solution as a function of time for both the true and assimilated systems.

5. sol_y_graphs
- Graphs in this directory display the y value solution as a function of time for both the true and assimilated systems.

6. sol_xy_graphs
- This directory contains graphs in which both the true and assimilated solutions are graphed separately

6a. sol_xy_graphs/animations
- This directory contains an animation of the graphs displayed in the directory "sol_xy_graphs".
- Due to the amount of time it takes to export one video, this application only renders the first matrix in the parent direcory 
  (I.e. sol_xy_graph_ ... _mtrx_0.jpg)

7. sol_xy_over_t_graphs
- Graphs in this directory display the true and assimilated x and y solutions displayed in two separate plots in one graph. 


Combined Graphs: 
- Graphs in this category provide a visual composite for a grouping of randomly generated matrices within a specific region of the trace-determinant.

1. tr_det_graph_ ... _cc.png
- A graph of the trace-determinant plane in which each randomly generated matrix is represented by a colored dot.
- The color of these dots are representative of each matrix's average relative error.  
- The color spectrum ranges from blue to red, in which completely blue dots have a low average relative error, and red dots do not. 

2. ev_graph_ ... _cc.png
- Follows the same dot-based coloring scheme of "tr_det_graph_ ... _cc.png". However, rather than using the trace-determinant plane, the dots here are instead graphed using their respective eigenvalues. 

3. pie_graph_ ... _cc.png
- Adheres the same coloring set-up of "tr_det_graph_ ... _cc.png".
 - However, instead of using dots, each partition of the graph represents the percentage of matrices that fall within a certain threshold range. 

4. bar_graph_ ... _cc.png
- Follows the same coloring set-up of "tr_det_graph_ ... _cc.png".
- However, each bar represents the frequency of threshold error category within each sample of matrices. 


Text Files:
1. ... _all_ matrices.txt
- A complete list of randomly generated matrices irregardless of their average relative error threshold. 

2. ... _bad_matrices.txt
- A list of "bad matrices" that fall below a certain average relative error threshold.

3. ... _good_matrices.txt
- A list of "good matrices" that are above a certain average relative error threshold.

--------------------------------------------------------------------------------------------------

OBSERVATIONS

When the nudging algorithm is applied to a given matrix, the average relative error for the parameter recovery process is below 10^−12. However, in two cases, the average relative error was between 10^−12 and 10^−8. To further analyze these occurrences, we will modify the code to pinpoint the exact matrices. 

ANTI-DIAGONAL
rde saddle 
rde sink
rde source

LEFT COLUMN
rde saddle 
rde source

MAIN DIAGONAL
ce sp_sink, 
rde saddle 
rde sink

RIGHT COLUMN 
ce sp source

--------------------------------------------------------------------------------------------------

Updated 2023.08.30 - 16:54
