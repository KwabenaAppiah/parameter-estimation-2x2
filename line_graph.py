#----- EXPORT TO IMAGE SEQUENCE ------------------#
# Standard library imports
import datetime
import time
import math
import os
from contextlib import ExitStack
import threading  # Added for timer
import subprocess  # For running external commands (Mov export)

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
# import imageio.v2 as imageio  # Import ImageIO v2 to avoid the Deprecation warning
import imageio
import shutil


class LineGraph:
    def __init__(self, *args):
        ev_type, pp_type, bounds, loop_limit, case_type = args

        # FOR COMP ONLY
        self._true_params = []
        self._has_non_optimal_matrices = False
        self._static_vars_dict = {
            'ev_type': ev_type, 'pp_type': pp_type, 'bounds': str(bounds), 'loop_limit': loop_limit, "case_type": case_type, "param_label_1": "", "param_label_2": ""

        }
        self.set_param_labels()
        # self.set_subplots_comp() 2023.09.09 hidden

    def set_param_labels(self):
        if self.get_static_vars_dict_elt("case_type") == "main_diagonal":
            self.set_static_vars_dict_elt("param_label_1", "$a_{11}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{22}$")

        elif self.get_static_vars_dict_elt("case_type") == "anti-diagonal":
            self.set_static_vars_dict_elt("param_label_1", "$a_{12}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{21}$")

        elif self.get_static_vars_dict_elt("case_type") == "left_column":
            self.set_static_vars_dict_elt("param_label_1", "$a_{11}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{21}$")

        elif self.get_static_vars_dict_elt("case_type") == "right_column":
            self.set_static_vars_dict_elt("param_label_1", "$a_{12}$")
            self.set_static_vars_dict_elt("param_label_2", "$a_{22}$")


    def organize_data(self, t, update_times, A, param_estimates, param_errors, avg_param_errors, curr_indx, non_optimal_threshold, S_lists, U_lists):
        U_lists_abs_err, U_lists_rel_err = U_lists
        # Init. static vars
        ev_type   = self.get_static_vars_dict_elt('ev_type')
        pp_type   = self.get_static_vars_dict_elt('pp_type')
        bounds    = self.get_static_vars_dict_elt('bounds')
        case_type = self.get_static_vars_dict_elt('case_type')
        loop_limit = self.get_static_vars_dict_elt('loop_limit')
        param_labels = (self.get_static_vars_dict_elt("param_label_1"), self.get_static_vars_dict_elt("param_label_2"))
        static_args = (ev_type, pp_type, bounds, loop_limit, param_labels, case_type)

        curr_indx = str(curr_indx)
        optimality_type = ""
        was_content_plotted = False

        # Parameter stuff
        param_1_estimates, param_2_estimates = param_estimates
        a11, a12, a21, a22  = A[0,0], A[0, 1], A[1,0], A[1,1]
        param_1_abs_err, param_2_abs_err, param_1_rel_err, param_2_rel_err = param_errors
        true_params = (a11, a12, a21, a22)
        true_params_str   = str(true_params)
        true_params_title = self.get_true_params_title(true_params)
        param_1 = param_2 = 0
        param_l_list, param_2_estimates = [], []
        avg_abs_param_err, avg_rel_param_err =  avg_param_errors

        linthresh_value = .000000000000001
        base_y = "y"
        base_x = "x"
        tilde_char = "\u0303"  # Unicode combining tilde character
        # Combine the base character and the combining tilde character
        yt = base_y + tilde_char
        xt = base_x + tilde_char
        assim_sol_chars = (xt, yt)


        if case_type == "main_diagonal":
            param_1, param_2 = a11, a22

        elif case_type == "anti-diagonal":
            param_1, param_2 = a12, a21

        elif case_type == "left_column":
            param_1, param_2 = a11, a21

        elif case_type == "right_column":
            param_1, param_2 = a12, a22


        error_type = ["absolute", "relative"]
        i = 0
        while i < len(error_type):

            if error_type[i] == "absolute":
                nth_avg_param_err = avg_abs_param_err[-1]
                U_list_err        = U_lists_abs_err
                param_1_err, param_2_err = param_1_abs_err, param_2_abs_err
                avg_param_error = avg_abs_param_err

            # else:
            elif error_type[i] == "relative":
                nth_avg_param_err = avg_rel_param_err[-1]
                U_list_err = U_lists_rel_err
                param_1_err, param_2_err = param_1_rel_err, param_2_rel_err
                avg_param_error = avg_rel_param_err

            if nth_avg_param_err >= non_optimal_threshold:
                optimality_type = "non-optimal"
                self.set_has_non_optimal_matrices_comp(True)
                was_content_plotted = True

            elif nth_avg_param_err < non_optimal_threshold:
                optimality_type = "optimal"
                was_content_plotted = True

            else:
                was_content_plotted = False

            # Error-based graphs
            self.display_avg_param_err(t, static_args, true_params_title, avg_param_error, curr_indx, linthresh_value, error_type[i])
            self.display_sep_param_err(t, static_args, param_1_err, param_2_err, true_params_title, curr_indx, linthresh_value,  error_type[i])
            self.display_sol_signal_err(t, static_args, U_list_err, true_params_title, curr_indx, linthresh_value,  error_type[i])
            self.display_sol_signal_err_split(t, static_args,  U_list_err, true_params_title, curr_indx, linthresh_value, error_type[i])
            self.write_to_file("mtrx_" + curr_indx + "|" + true_params_str, case_type, optimality_type, error_type[i])
            i = i + 1

        # Non-error-based graphs
        # self.write_to_file("mtrx_" + curr_indx + "|" + true_params_str, case_type, "all", "N/A")
        self.write_to_file("mtrx_" + curr_indx + "|" + true_params_str, case_type)
        self.display_sol_xy_over_t(t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars)
        self.display_sol_x(t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars)
        self.display_sol_y(t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars)
        self.display_sol_xy(static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars, False)

        return was_content_plotted

    def write_to_file(self, line_str, case_type, optimality_type = "all", error_type = "N/A"):
        if error_type == "absolute":
            err_label = "abs"
        else:
            err_label = "rel"

        output = line_str + "\n" + "\n"
        ev_and_pp_type = self.get_static_vars_dict_elt("ev_type") + "_" + self.get_static_vars_dict_elt("pp_type")
        if error_type == "N/A":
            subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/text_files/"
            filename = ev_and_pp_type + "_" + optimality_type + "_matrices.txt"
        else:
            subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/text_files/" + error_type + "_error/"
            filename = ev_and_pp_type + "_" + err_label  + "_err_" + optimality_type + "_matrices.txt"



        os.makedirs(subdir, exist_ok = True)
        f = open(subdir + filename, "a")
        f.write(output)
        f.close()

    def format_fl_vals(self, val):
        return "{:.5f}".format(val)


    def display_avg_param_err(self, t, static_args, true_params_title, avg_param_err, curr_indx, linthresh_value, error_type):
        err_label_1 = err_label_2 = ""

        if(error_type == "absolute"):
            err_label = "abs"
            line_color = "tomato"

        else:
            err_label = "rel"
            line_color = "royalblue"

        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        param_label_1, param_label_2 = param_labels

        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.plot(t, avg_param_err, color = line_color, label = param_label_1 + " & " + param_label_2 + ": \n" + "Avg. " + err_label + " err" )
        ax.set_xlabel("Time")
        ax.set_ylabel("Avg. " + error_type + " parameter error - " + param_label_1 + " & " +  param_label_2)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh = linthresh_value)

        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type.upper()  + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        ax.set_title(label = title, pad = 20, fontsize = 15)

        # Output image files
        # output_type = "indv"
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/avg_param_err_graphs/" + error_type + "_error/"
        filename = "avg_" + err_label + "_param_err_graph_" + ev_and_pp_type
        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".jpg", dpi = 300)
        plt.close('all')


    def display_sep_param_err(self, t, static_args, param_1_err, param_2_err, true_params_title, curr_indx, linthresh_value, error_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        param_label_1, param_label_2 = param_labels
        err_label_1 = err_label_2 = ""
        if(error_type == "absolute"):
            err_label =  "abs"
            line_color_1, line_color_2 = "tomato", "yellowgreen"

        else:
            err_label = "rel"
            line_color_1, line_color_2 = "royalblue", "orange"

        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.set_xlabel("Time")
        ax.set_ylabel(error_type.capitalize() + " parameter error - " + param_label_1 + " & " + param_label_2)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh = linthresh_value)
        ax.plot(t, param_1_err, color = line_color_1, label = param_label_1 + " " + err_label + ". err." )
        ax.plot(t, param_2_err, color = line_color_2, label = param_label_2 + " " + err_label + ". err." )

        # For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        ax.set_title(label = title, pad = 20, fontsize = 15)

        # Output img files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() +  "/sep_param_err_graphs/" + error_type +  "_error/"

        filename = "sep_"+ err_label + "_param_err_graph_" + ev_and_pp_type
        ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".jpg", dpi = 300)
        plt.close('all')


    def display_sol_xy(self, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars, is_animation_on):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, (ax_1, ax_2) = plt.subplots(1, 2)
        fig.set_size_inches(30, 15)
        label_font_size = 18

        line_color_1, line_color_2 = "yellowgreen", "orange"
        xt, yt = assim_sol_chars
        ax_1.set_xlabel("x", fontsize = label_font_size)
        ax_1.set_ylabel("y", fontsize = label_font_size)
        ax_1.plot(x_estimates, y_estimates, label = "True sol.", color = line_color_1)
        ax_1.set_xscale("symlog", linthresh = linthresh_value)
        ax_1.set_yscale("symlog", linthresh = linthresh_value)
        ax_1.legend(loc = "upper right")

        ax_2.set_xlabel(xt, fontsize = label_font_size)
        ax_2.set_ylabel(yt, fontsize = label_font_size)
        ax_2.plot(xt_estimates, yt_estimates, label = "Est. sol.", color = line_color_2)
        ax_2.set_xscale("symlog", linthresh = linthresh_value)
        ax_2.set_yscale("symlog", linthresh = linthresh_value)
        ax_2.legend(loc = "upper right")

        # For file title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = label_font_size)
        title_specs = (title, label_font_size)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type

        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/sol_xy_graphs/"
        filename = "sol_xy_graph" + "_" + ev_and_pp_type + "_mtrx_" + curr_indx

        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + ".jpg", dpi = 300)

        # if(is_animation_on == True and int(curr_indx) < 1 or optimality_type == "non-optimal" and is_animation_on == True): #Only output an animation of the first matrix or a non-optimal matrix
        if(is_animation_on == True and int(curr_indx) < 1): #Only output an animation of the first matrix or a non-optimal matrix
            animation_subdir = subdir + "animations/"
            os.makedirs(animation_subdir, exist_ok = True)
            self.animation_output_1x2(x_estimates, y_estimates, xt_estimates, yt_estimates, linthresh_value, animation_subdir + filename, assim_sol_chars, title_specs)
        plt.close('all')




    # ----- EXPORT TO IMAGE SEQUENCE ------------------#
    def animation_output_1x2(self, x1, y1, x2, y2, linthresh_value, video_filename, assim_sol_chars, title_specs):
        is_exporting_to_mp4 = True
        is_exporting_to_gif = False
        xt, yt = assim_sol_chars
        ev_type   = self.get_static_vars_dict_elt('ev_type')

        # Initialize the figure and subplots
        fig, axes = plt.subplots(1, 2, figsize = (30, 15))

        # Calculate the number of frames for the animation
        num_frames = len(x1)

        # Create a directory to store the frames
        frames_dir = "temp_animation_frames"
        line_color_1, line_color_2 = "yellowgreen", "orange"


        # Delete the directory if it already exists and recreate it
        if os.path.exists(frames_dir):
            shutil.rmtree(frames_dir)
        os.makedirs(frames_dir)

        # Generate frames for the animation

        j = 0
        if(ev_type =="ce"):
            mod_val = 20 #20 for the standard loop, # Set to 400 for debugging purposes, as it will provide for a shorter loop
        else:
            mod_val = 100

        for i in range(num_frames):
            # Set up plots
            if i % mod_val == 0:
                axes[0].set_title('True sol.')
                axes[0].set_xlabel('x')
                axes[0].set_ylabel('y')
                axes[0].set_xscale('symlog', linthresh = linthresh_value)
                axes[0].set_yscale('symlog', linthresh = linthresh_value)
                axes[0].plot(x1[:i + 1], y1[:i + 1], color = line_color_1, label = "True sol." )
                axes[0].set_xlim(np.min(x1), np.max(x1))
                axes[0].set_ylim(np.min(y1), np.max(y1))

                axes[1].set_title('Est. sol.')
                axes[1].set_xlabel(xt)
                axes[1].set_ylabel(yt)
                axes[1].set_xscale('symlog', linthresh = linthresh_value)
                axes[1].set_yscale('symlog', linthresh = linthresh_value)
                axes[1].plot(x2[:i + 1], y2[:i + 1], color = line_color_2, label = "Est. sol.")
                axes[1].set_xlim(np.min(x2), np.max(x2))
                axes[1].set_ylim(np.min(y2), np.max(y2))

                if(i == 0): #Note: without this line of code, the axes[0].legend and axex[1].legend would print out each iteration
                    axes[0].legend(loc = "upper right")
                    axes[1].legend(loc = "upper right")
                    title, label_font_size = title_specs
                    fig.suptitle(title, fontsize = label_font_size)

                # Save the current plot as an image
                frame_file = os.path.join(frames_dir, f"frame_{j:04d}.jpg")
                plt.savefig(frame_file, format = 'jpg')
                print("Export frame:", j + 1, '/', int(num_frames / mod_val), end='\r')
                j += 1

        # Close the current figure
        plt.close()
        print("")

        # Loop through the range from 0 to j (exclusive)
        frames = []
        for i in range(j):
            # Construct the filename using string formatting to pad the number with zeros
            filename = os.path.join(frames_dir, f"frame_{i:04d}.jpg")

            # Read the image using imageio.v2.imread and append it to the frames list
            frames.append(imageio.imread(filename))
            print("Composited frames:", i + 1, '/', j, end = '\r')

        if(is_exporting_to_mp4 == True):
            self.export_images_for_mp4(frames)
            self.convert_images_to_mp4(frames_dir, video_filename)

        elif(is_exporting_to_gif == True):
            self.export_to_gif(frames, frames_dir, video_filename)


 #********** FOR EXPORTING MP4s ***********************************/
    def export_images_for_mp4(self, frames):
        frames_dir = "temp_animation_frames"
        if os.path.exists(frames_dir):
            shutil.rmtree(frames_dir)
        os.makedirs(frames_dir)

        for i, frame in enumerate(frames):
            frame_file = os.path.join(frames_dir, f"frame_{i:04d}.jpg")
            imageio.imwrite(frame_file, frame)
            print("Exported frames:", i + 1, '/', len(frames), end ='\r')

    def convert_images_to_mp4(self, frames_dir, video_filename):
        mp4_filename = video_filename + '.mp4'
        cmd = ['ffmpeg', '-framerate', '10', '-i', os.path.join(frames_dir, 'frame_%04d.jpg'), '-c:v', 'libx264', '-r', '30', mp4_filename]
        subprocess.run(cmd)

        if os.path.exists(mp4_filename):
            print("\n" + "'" + mp4_filename + "'", "now exists.")
            shutil.rmtree(frames_dir) #Remove frame directory


  #********** FOR EXPORTING GIFS ***********************************/
    def export_to_gif(self, frames, frames_dir, gif_filename):
        print("")
        #### Timer
        stop_event = threading.Event()
        timer_thread = threading.Thread(target = self.init_timer, args = (stop_event,)) # With trailing comma since it's a tupple
        timer_thread.start()
        ####
        gif_filename = gif_filename + '.gif'
        imageio.mimsave(gif_filename, frames, duration = 0.001)  # 10 / 10000 = 0.001 a duration of 10 seconds
        if os.path.exists(gif_filename):
            print("\n" + "'" + gif_filename + "'", "now exists.")
            self.stop_timer(stop_event)
            timer_thread.join()

            # Remove the temporary frames directory
            shutil.rmtree(frames_dir)

        else:
            self.init_timer(stop_event)


    def format_time(self, seconds):
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)

        time_str = ""
        if days:
            time_str += f"{days} days, "

        elif hours:
            time_str += f"{hours} hours, "

        elif minutes:
            time_str += f"{minutes} minutes, "

        time_str += f"{seconds:.2f} seconds"

        return time_str

    def init_timer(self, stop_event):
        start_time = time.time()
        while not stop_event.is_set():
            elapsed_time = time.time() - start_time
            time_str = self.format_time(elapsed_time)
            print(f"\rFinal Rendering Time: {time_str}", end = "")
            time.sleep(0.1)  # Adjust the sleep time to control the timer's granularity

    def stop_timer(self, stop_event):
        stop_event.set()

    #********** For Exporting gifs - END ***********************************/
    def display_sol_xy_over_t(self, t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars):
    # def display_sol_xy_over_t(self, t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars, optimality_type):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, (ax_1, ax_2) = plt.subplots(2, 1)
        fig.set_size_inches(16, 12)
        xt, yt = assim_sol_chars
        line_color_1, line_color_2 = "yellowgreen", "orange"

        # For x and xt
        ax_1.set_ylabel("x and " + xt)
        ax_1.set_xlabel("Time")
        ax_1.plot(t, x_estimates, label = "x", color = line_color_1)
        ax_1.plot(t, xt_estimates, label = xt, color = line_color_2)
        ax_1.legend(loc = "best")
        ax_1.set_xscale("log")
        ax_1.set_yscale("symlog", linthresh = linthresh_value)

        # For y and yt
        ax_2.set_ylabel("y and " + yt)
        ax_2.set_xlabel("Time")
        ax_2.plot(t, y_estimates, label = "y", color = line_color_1)
        ax_2.plot(t, yt_estimates, label = yt, color = line_color_2)
        ax_2.legend(loc = "best")
        ax_2.set_xscale("log")
        ax_2.set_yscale("symlog", linthresh = linthresh_value)


        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/sol_xy_over_t_graphs/"
        filename = "sol_xy_over_t_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".jpg", dpi = 300)
        plt.close('all')


    def display_sol_x(self, t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        line_color_1, line_color_2 = "yellowgreen", "orange"

        # For x and xt
        xt, yt = assim_sol_chars
        ax.set_ylabel("x and " + xt)
        ax.set_xlabel("Time")
        ax.plot(t, x_estimates, label = 'x', color = line_color_1)
        ax.plot(t, xt_estimates, label = xt, color = line_color_2)
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh = linthresh_value)
        ax.legend(loc = "best")
        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_"+ self.get_date_str() + "/sol_x_graphs/"
        filename = "sol_x_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".jpg", dpi = 300)
        plt.close('all')


    def display_sol_y(self, t, static_args, S_lists, true_params_title, curr_indx, linthresh_value, assim_sol_chars):
        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_estimates, y_estimates, xt_estimates, yt_estimates = S_lists
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        line_color_1, line_color_2 = "yellowgreen", "orange"

        # For y and yt
        xt, yt = assim_sol_chars
        ax.set_ylabel("y and " + yt)
        ax.set_xlabel("Time")
        ax.plot(t, y_estimates, label = 'y', color = line_color_1)
        ax.plot(t, yt_estimates, label = yt, color = line_color_2)
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh = linthresh_value)
        ax.legend(loc = "best")

        title = ev_type.upper() + " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/sol_y_graphs/"
        filename = "sol_y_graph" + "_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".jpg", dpi = 300)
        plt.close('all')


    def display_sol_signal_err(self, t, static_args, U_lists, true_params_title, curr_indx, linthresh_value, error_type):
        if(error_type == "absolute"):
            err_label =  "abs"
            line_color_1, line_color_2 = "tomato", "yellowgreen"
        else:
            err_label = "rel"
            line_color_1, line_color_2 = "royalblue", "orange"

        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_sig_err_estimates, y_sig_err_estimates = U_lists
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 8)
        ax.set_ylabel(error_type.capitalize() +  " signal error")
        ax.set_xlabel("Time")
        ax.plot(t, x_sig_err_estimates, label = "x " + err_label + ". signal err.", color = line_color_1)
        ax.plot(t, y_sig_err_estimates, label = "y " + err_label + ". signal err.", color = line_color_2)
        ax.legend(loc = "best")
        ax.set_xscale("log")
        ax.set_yscale("log")

        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/sol_signal_err_graphs/" + error_type + "_error/"
        filename =  "sol_" + err_label + "_signal_err_graph_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".jpg", dpi = 300)
        plt.close('all')



    def display_sol_signal_err_split(self, t, static_args, U_lists, true_params_title, curr_indx, linthresh_value, error_type):
        if(error_type == "absolute"):
            err_label =  "abs"
            line_color_1, line_color_2 = "tomato", "yellowgreen"
        else:
            err_label = "rel"
            line_color_1, line_color_2 = "royalblue", "orange"

        ev_type, pp_type, bounds, loop_limit,  param_labels, case_type = static_args
        x_sig_err_estimates, y_sig_err_estimates = U_lists
        fig, (ax_1, ax_2) = plt.subplots(2, 1)
        fig.set_size_inches(16, 12)

        ax_1.set_ylabel( error_type.capitalize() + " signal error")
        ax_1.set_xlabel("Time")
        ax_1.plot(t, x_sig_err_estimates, label = "x " + err_label  + ". signal err.", color = line_color_1 )

        ax_1.legend(loc = "best")
        ax_1.set_xscale("log")
        ax_1.set_yscale("log")

        ax_2.set_ylabel(error_type.capitalize() + " signal error")
        ax_2.set_xlabel("Time")
        ax_2.plot(t, y_sig_err_estimates, label = "y " + err_label  + ". signal err.", color = line_color_2)
        ax_2.legend(loc = "best")
        ax_2.set_xscale("log")
        ax_2.set_yscale("log")

        #For title
        title = ev_type.upper()+ " | " + pp_type.upper() + " | " + case_type.upper() + " | BNDS " + bounds + " | MTRX " + curr_indx + " : " + true_params_title
        fig.suptitle(title, fontsize = 15)

        # Output IMG files
        ev_and_pp_type = ev_type + "_" + pp_type
        subdir = "../output/" + ev_and_pp_type + "_" + case_type + "_" + self.get_date_str() + "/sol_signal_err_split_graphs/" + error_type + "_error/"
        filename = "sol_" + err_label + "_signal_err_split_graph_" + ev_and_pp_type
        os.makedirs(subdir, exist_ok = True)
        fig.savefig(subdir + filename + "_mtrx_" + curr_indx + ".jpg", dpi = 300)
        plt.close('all')


#-------------- GETTERS ---------------------------------------------


    def get_date_str(self):
        raw_date = datetime.datetime.now()
        t = time.localtime()
        date_str = str(raw_date.year) + "."+ str(time.strftime("%m"))+ "." + str(time.strftime("%d"))
        return date_str

    def get_subplots_comp(self):
        return self._fig, self._ax

    def get_true_params_comp(self):
        return self._true_params

    def has_non_optimal_matrices_comp(self):
        return self._has_non_optimal_matrices

    def get_static_vars_dict(self):
        return self._static_vars_dict

    def get_static_vars_dict_elt(self, key):
        return self._static_vars_dict[key]


#-------------- SETTERS -------------------------------------
    def set_has_non_optimal_matrices_comp(self, val):
        self._has_non_optimal_matrices = val

    def set_ev_type(self, val):
        self._ev_type = val

    def set_pp_type(self, val):
        self._pp_type = val

    def set_static_vars_dict_elt(self, key, value):
        self._static_vars_dict[key] = value


    def get_true_params_title(self, true_params):
        a11, a12, a21, a22 = true_params
        a11_str, a12_str = self.format_fl_vals(a11), self.format_fl_vals(a12)
        a21_str, a22_str = self.format_fl_vals(a21), self.format_fl_vals(a22)
        return "[" + a11_str + ", " + a12_str + ", " + a21_str + ", " + a22_str + "]"

    # def animate_loading(self):
    #     loading_text = "Loading"
    #     while True:
    #         for i in range(4):
    #             sys.stdout.write(f"\r{loading_text}{'.' * i}   ")
    #             sys.stdout.flush()
    #             time.sleep(0.5)
