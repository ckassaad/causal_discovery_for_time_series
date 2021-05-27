if __name__ == "__main__":
    measure_name = "F-Score"  # Precision Recall F-Score

    structure_list = ['dsr_fork', 'dsr_v_structure', 'dsr_mediator', 'dsr_diamond']
    structure_hidden = []
    n_samples_list = [125, 250, 500, 1000]

    path_output = "../experiments/performance_average/latex_format/figures" + "/" + measure_name + \
                  "_dsr_summarySelf_summaryExternal_latex_figure_.txt"
    f_o = open(path_output, "w")
    f_o.write("\\begin{figure}%[ht!]\n\t\centering\n")

    option_by_method = {'PCTMIwindow=auto': 'blue,smooth,mark=*',
                       #  'NBCB': 'red,smooth,mark=*',
                       # 'NBCB_pw' : 'red,smooth,mark=o',
                       #  'PWNBCBk': 'purple,smooth,mark=*',
                        'GrangerPW': 'black,smooth,mark=o',
                        'GrangerMV': 'black,smooth,mark=*',
                        # 'GrangerK': 'black,smooth,mark=x',
                        'TCDF': 'gray,smooth,mark=*',
                        'PCMCICMIknn': 'brown, smooth,mark=*',
                        'PCMCIParCorr': 'brown, smooth,mark=o',
                        # 'oCSE': 'green,smooth,mark=*',
                        # 'tsFCI': 'pink,smooth,mark=*',
                        # 'FCITMI': 'cyan,smooth,mark=*',
                        # 'VarLiNGAM': 'yellow,smooth,mark=*',
                        'TiMINO': 'orange,smooth,mark=*',
                        'Dynotears': 'pink,smooth,mark=*'}
    add_subtract_from_sample_size_by_method = {'PCTMIwindow=auto': -10,
                                               #  'NBCB': -9,
                                               # 'NBCB_pw' : -8,
                                               # 'PWNBCBk' : -7,
                                               'GrangerPW': -6,
                                                'GrangerMV': -5,
                                                # 'GrangerK': -4,
                                                'TCDF': -3,
                                                'PCMCICMIknn': -2,
                                                'PCMCIParCorr': -1,
                                                # 'oCSE': 1,
                                                # 'tsFCI': 2,
                                                # 'FCITMI': 3,
                                                # 'VarLiNGAM': 4,
                                                'TiMINO': 5,
                                                'Dynotears': 6}
    for structure in structure_list:
        measure = ""
        for type_of_causes in ["summary_other", "summary_other", "summary_self"]:
            path_input = "../experiments/performance_average/" + type_of_causes + "_performance_average"
            if type_of_causes == "summary_other":
                # all_methods_list = ['PCTMIwindow=auto', 'NBCB', 'NBCB_pw', 'PWNBCBk', 'GrangerPW', 'GrangerMV', 'TCDF',
                #                     'PCMCICMIknn', 'PCMCIParCorr', 'oCSE', 'VarLiNGAM', 'TiMINO']
                all_methods_list = ['PCTMIwindow=auto', 'GrangerPW', 'GrangerMV', 'TCDF',
                                    'PCMCICMIknn', 'PCMCIParCorr', 'TiMINO', 'Dynotears']
                methods_treats_hidden = []
                measure_type1 = "Adjacent"#
                measure_type2 = "Oriented"#
                if measure == measure_name + " " + measure_type1:#
                    measure = measure_name + " " + measure_type2#
                else:#
                    measure = measure_name + " " + measure_type1#
                    # measure_type = "Oriented"
                    # measure = measure_name + " " + measure_type
                    legend = "\t\t\\legend{"
                    addplot_for_legend = ""
                    for method in all_methods_list:
                        legend = legend + "{" + method + "},"
                        addplot_for_legend = addplot_for_legend + "\t\t\\addplot[" + option_by_method[
                            method] + ",] coordinates {(0, -1)};\n"
                    legend = legend[:-1]
                    legend = legend + "}\n"
            elif type_of_causes == "summary_self":
                # all_methods_list = ['PCTMIwindow=auto', 'NBCB', 'NBCB_pw', 'PWNBCBk', 'oCSE', 'TCDF', 'PCMCICMIknn', 'PCMCIParCorr', 'VarLiNGAM']
                all_methods_list = ['PCTMIwindow=auto', 'TCDF', 'PCMCICMIknn', 'PCMCIParCorr', 'Dynotears']
                methods_treats_hidden = []
                measure = measure_name
            else:
                all_methods_list = None
                methods_treats_hidden = None

            if structure in structure_hidden:
                methods_list = methods_treats_hidden
            else:
                methods_list = all_methods_list
            if structure == structure_list[0]:
                if type_of_causes == "summary_other":
                    if "Oriented" in measure:
                        title = "\n\t\ttitle = {Summary: External Causation},"
                    else:
                        title = "\n\t\ttitle = {Summary: External Adjacency},"
                else:
                    title = "\n\t\ttitle = {Summary: Self Causation},"
            else:
                title = ""
            if (type_of_causes == "summary_other") and (measure == measure_name + " " + "Adjacent"):
                subfigure_begin = "\t\\resizebox {0.32\\textwidth} {!} {\n\t\\begin{subfigure}[b]{0.5\\textwidth}\n\t\t\\begin{tikzpicture}[font=\small]\n\t\t\\renewcommand{\\axisdefaulttryminticks}{4}\n\t\t\pgfplotsset{every major grid/.append style={densely dashed}}\n\t\t\pgfplotsset{every axis legend/.append style={cells={anchor=west},fill=white, at={(0.02,0.98)}, anchor=north west}}\n\t\t\\begin{axis}["+title+"\n\t\txmin = 100,\n\t\txmax = 1100,\n\t\txmode=log,\n\t\tlog ticks with fixed point,\n\t\tymin=0,\n\t\tymax=1,\n\t\tgrid=minor,\n\t\tscaled ticks=true,\n\t\tylabel = {" + structure.replace("_", "-") + "},\n\t\theight = 4.5cm,\n\t\twidth=7cm,\n\t\tlegend style={nodes={scale=0.65, transform shape}}\n\t\t]\n"
            else:
                subfigure_begin = "\t\\resizebox {0.32\\textwidth} {!} {\n\t\\begin{subfigure}[b]{0.5\\textwidth}\n\t\t\\begin{tikzpicture}[font=\small]\n\t\t\\renewcommand{\\axisdefaulttryminticks}{4}\n\t\t\pgfplotsset{every major grid/.append style={densely dashed}}\n\t\t\pgfplotsset{every axis legend/.append style={cells={anchor=west},fill=white, at={(0.02,0.98)}, anchor=north west}}\n\t\t\\begin{axis}["+title+"\n\t\txmin = 100,\n\t\txmax = 1100,\n\t\txmode=log,\n\t\tlog ticks with fixed point,\n\t\tymin=0,\n\t\tymax=1,\n\t\tgrid=minor,\n\t\tscaled ticks=true,\n\t\theight = 4.5cm,\n\t\twidth=7cm,\n\t\tlegend style={nodes={scale=0.65, transform shape}}\n\t\t]\n"
            subfigure_end = "\t\t\end{axis}\n\t\t\end{tikzpicture}\n\t\t\caption{}\n\t\t\label{}\n\t\end{subfigure}\n}\n\hfill\n"
            f_o.write(subfigure_begin)
            for method in methods_list:
                add_subtract = add_subtract_from_sample_size_by_method[method]
                addplot_begin = "\t\t\\addplot["+option_by_method[method]+", error bars/.cd, y dir=both,y explicit] plot coordinates{\n"
                addplot_end = "\t\t};\n"
                f_o.write(addplot_begin)
                for n_samples in n_samples_list:
                    add_subtract_multiply = int(n_samples>1000)*10
                    add_subtract_multiply = 1 if add_subtract_multiply == 0 else add_subtract_multiply
                    files_input_name = method + "_" + structure + "_" + str(n_samples)
                    try:
                        f_i = open(str(path_input) + "/" + str(files_input_name), "r")
                        for line_i in f_i:
                            if measure in line_i:
                                nextLine = next(f_i)
                                # nextLine = nextLine.replace("+-", "\pm")
                                print(method, structure, n_samples, type_of_causes)
                                print(nextLine)
                                nextLine = nextLine.replace("\n", "")
                                nextline_processed = nextLine.split(" ")
                                measure_mean = "{:.2f}".format(round(float(nextline_processed[0]), 2))
                                measure_std = "{:.2f}".format(round(float(nextline_processed[2]), 2))
                                result_line = "\t\t\t("+str(n_samples+add_subtract*add_subtract_multiply)+", "+measure_mean+") +- ("+measure_std+", "+measure_std+")\n"
                                print(result_line)
                                f_o.write(result_line)
                                break
                        f_i.close()
                    except FileNotFoundError:
                        print(files_input_name)
                        # f_o.write("f\n")
                f_o.write(addplot_end)
            f_o.write(subfigure_end)
    legend = legend.replace("PCTMIwindow=auto", "PCTMI")
    # legend = legend.replace("NBCB_pw", "pwNBCB")
    # legend = legend.replace("PWNBCBk", "pwNBCBk")
    legend = legend.replace("Granger", "GC")
    legend = legend.replace("PCMCICMIknn", "PCMCIMI")
    legend = legend.replace("PCMCIParCorr", "PCMCIPC")
    legend = legend.replace("TiMINO", "TiMINo")
    outside_legend = "\\begin{tikzpicture}[baseline=(leg.center)]\\begin{axis}[width=0.2\\textwidth, height=0.2 \\textwidth,hide axis,xmin=0, xmax=10, ymin=0, ymax=0.0,legend columns=6,legend style={name=leg,draw=white!15!black,legend cell align=left,at={(0.5,0.5)},font=\\tiny}]" + addplot_for_legend + legend + "\\end{axis}\\end{tikzpicture}\n"
    f_o.write(outside_legend)
    f_o.write("\\caption"
            "[Detailed perfomances]{Adjacency (F1), external causation ($\overrightarrow{\text{F1}} $) and self causation ($\mathring{\text{F1}}$) in the summary causation graph for all the methods on 4 simulated datasets (mean $\pm$ standard deviation). Results  are computed for various time grid sizes, from 125 to 1000. A log-scale is used in abscissa.}"
              "\n\label{fig:simu:details_dsr}\n\end{figure}")
    f_o.close()

