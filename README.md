# parameter-estimation

BACKGROUND

Based on the findings of Elizabeth Carlson, Joshua Hudson and Adam Larios (https://epubs.siam.org/doi/10.1137/19M1248583),
this project utilizes a nudging algorithm for the purpose of parameter recovery. As of now, this script only works for two-dimensional linear systems.

In addition, this project is a modified version of code originally developed by Nathan Taylor
under the guidance of Dr. Vincent R Martinez at the Department of Mathematics and Statistics at Hunter College, City University of New York
(http://math.hunter.cuny.edu/vmartine/).

For more information on how to properly setup an environment to execute this code
and its predecessor, please check here: https://github.com/taylorbn/linear_nudging.


--------------------------------------------------------------------------------------------------


HOW DO I RUN THIS SCRIPT?

SAMPLE COMMAND LINE PROMPTS:
- python parameter_estimation.py rde saddle 100 10 5 30
- python parameter_estimation.py ce sp_sink 90 10 7 25 

CONTEXT: 
- python parameter_estimation.py EV_TYPE PP_TYPE MU_VAL RLX_TIME BND_VAL LOOP_LIM

EV_TYPE (eigenvalue type) and PP_TYPE (phase portrait type):


1. rde = real distinct eigenvalues
- saddle
- source
- sink


2. re = repeated eigenvalues
- sink
- source

Note: Currently there’s an issue with the trace-determinant graph for these two areas.


3. ce =  complex eigenvalues
- sp_sink = Spiral sink, 
- sp_source = spiral source
- center


MU_VAL (µ)
- The nudging parameter

RLX_TIME: 
- Relaxation time

BND_VAL: 
- Bound Val, or the boundaries of which the randomly generated matrix values are pulled from.
Thus, the sample line above would cover a boundary from -5 to 5, inclusive.
Note: When the BND_VAL is increased, the algorithm starts to perform poorly.
Thus, as of now, a bound value of 5 is recommended for testing purposes.  


LOOP_LIM: 
- The total number of randomly generated matrices of which the nudging algorithm will be applied to.


--------------------------------------------------------------------------------------------------


WHAT DOES THIS SCRIPT OUTPUT?

This script outputs four types of documents into a folder called output (which is generated dynamically if non-existent).
Please note, the filenames may change depending on the eigenvalue and phase portrait of the matrices.

1. tr_det_graph.png

- A graph of the trace-determinant plane with each randomly generated matrix represented by a colored dot.
Moreover, the color of these dots are representative of each matrix’s average relative error. Currently, this is the only
graphic that visualizes all of the matrices generated.

2. bad_matricies.txt

- A list of “bad matricies” above a certain average relative error threshold.

3. line_graph_indv.png

- An individual line graph for each matrix with a certain average relative error threshold.

4. line_graph_comp.png

- A composite line graph that contains all of the bad matrices produced. Moreover, the content
found here should mirror that of the content found in bad_matricies.txt.


Updated 2023.03.13 - 04:18
