if __name__ == "__main__":
    measure_name = "Recall"  # Precision Recall F-Score

    structure_list = ['fork', 'v_structure', 'mediator', 'diamond', '7ts2h']
    structure_hidden = ['7ts2h']
    n_samples_list = [125, 250, 500, 1000, 2000, 4000]

    path_output = "../experiments/performance_average/latex_format/figures" + "/" + measure_name + \
                  "_survey_summary_details_latex_figure.txt"
    f_o = open(path_output, "w")
    f_o.write("\\begin{figure}%[ht!]\n\t\centering\n")

    option_by_method = {'GrangerPW': 'black,smooth,mark=o',
                        'GrangerMV': 'black,smooth,mark=*',
                        'GrangerK': 'black,smooth,mark=x',
                        'TCDF': 'gray,smooth,mark=*',
                        'PCMCICMIknn': 'brown, smooth,mark=*',
                        'PCMCIParCorr': 'brown, smooth,mark=o',
                        'oCSE': 'green,smooth,mark=*',
                        'tsFCI': 'pink,smooth,mark=*',
                        'VarLiNGAM': 'yellow,smooth,mark=*',
                        'TiMINO': 'orange,smooth,mark=*',
                        'Dynotears': 'purple,smooth,mark=*'}
    add_subtract_from_sample_size_by_method = {'GrangerPW': -6,
                                                'GrangerMV': -5,
                                                'GrangerK': -4,
                                                'TCDF': -3,
                                                'PCMCICMIknn': -2,
                                                'PCMCIParCorr': -1,
                                                'oCSE': 1,
                                                'tsFCI': 2,
                                                'VarLiNGAM': 4,
                                                'TiMINO': 5,
                                               'Dynotears': 6}
    for structure in structure_list:
        measure = ""
        for type_of_causes in ["summary_other", "summary_self"]:
            path_input = "../experiments/performance_average/" + type_of_causes + "_performance_average"
            if type_of_causes == "summary_other":
                all_methods_list = ['GrangerPW', 'GrangerMV', 'TCDF', 'PCMCICMIknn', 'PCMCIParCorr', 'oCSE', 'tsFCI', 'VarLiNGAM', 'TiMINO', 'Dynotears']
                methods_treats_hidden = ['TCDF', 'tsFCI']
                # measure_type1 = "Adjacent"
                measure_type = "Oriented"
                measure = measure_name + " " + measure_type
                legend = "\t\t\\legend{"
                addplot_for_legend = ""
                for method in all_methods_list:
                    legend = legend + "{" + method + "},"
                    addplot_for_legend = addplot_for_legend + "\t\t\\addplot[" + option_by_method[
                        method] + ",] coordinates {(0, -1)};\n"
                legend = legend[:-1]
                legend = legend + "}\n"
            elif type_of_causes == "summary_self":
                all_methods_list = ['TCDF', 'PCMCICMIknn', 'PCMCIParCorr', 'oCSE', 'tsFCI', 'VarLiNGAM', 'Dynotears']
                methods_treats_hidden = ['TCDF', 'tsFCI']
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
                    title = "\n\t\ttitle = {Summary: without self causes},"
                else:
                    title = "\n\t\ttitle = {Summary: only self causes},"
            else:
                title = ""
            if type_of_causes == "summary_other":
                subfigure_begin = "\t\\resizebox {0.45\\textwidth} {!} {\n\t\\begin{subfigure}[b]{0.5\\textwidth}\n\t\t\\begin{tikzpicture}[font=\small]\n\t\t\\renewcommand{\\axisdefaulttryminticks}{4}\n\t\t\pgfplotsset{every major grid/.append style={densely dashed}}\n\t\t\pgfplotsset{every axis legend/.append style={cells={anchor=west},fill=white, at={(0.02,0.98)}, anchor=north west}}\n\t\t\\begin{axis}["+title+"\n\t\txmin = 0,\n\t\txmax = 4500,\n\t\txmode=log,\n\t\tlog ticks with fixed point,\n\t\tymin=0,\n\t\tymax=1,\n\t\tgrid=minor,\n\t\tscaled ticks=true,\n\t\tylabel = {" + structure.replace("_", "-") + "},\n\t\theight = 4.5cm,\n\t\twidth=8cm,\n\t\tlegend style={nodes={scale=0.65, transform shape}}\n\t\t]\n"
            else:
                subfigure_begin = "\t\\resizebox {0.45\\textwidth} {!} {\n\t\\begin{subfigure}[b]{0.5\\textwidth}\n\t\t\\begin{tikzpicture}[font=\small]\n\t\t\\renewcommand{\\axisdefaulttryminticks}{4}\n\t\t\pgfplotsset{every major grid/.append style={densely dashed}}\n\t\t\pgfplotsset{every axis legend/.append style={cells={anchor=west},fill=white, at={(0.02,0.98)}, anchor=north west}}\n\t\t\\begin{axis}["+title+"\n\t\txmin = 0,\n\t\txmax = 4500,\n\t\txmode=log,\n\t\tlog ticks with fixed point,\n\t\tymin=0,\n\t\tymax=1,\n\t\tgrid=minor,\n\t\tscaled ticks=true,\n\t\theight = 4.5cm,\n\t\twidth=8cm,\n\t\tlegend style={nodes={scale=0.65, transform shape}}\n\t\t]\n"
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
    legend = legend.replace("Granger", "GC")
    legend = legend.replace("PCMCICMIknn", "PCMCIMI")
    legend = legend.replace("PCMCIParCorr", "PCMCIPC")
    outside_legend = "\\begin{tikzpicture}[baseline=(leg.center)]\\begin{axis}[width=0.2\\textwidth, height=0.2 \\textwidth,hide axis,xmin=0, xmax=10, ymin=0, ymax=0.0,legend columns=6,legend style={name=leg,draw=white!15!black,legend cell align=left,at={(0.5,0.5)},font=\\tiny}]" + addplot_for_legend + legend + "\\end{axis}\\end{tikzpicture}\n"
    f_o.write(outside_legend)
    f_o.write("\\caption"
              "[Global performance]{Global performance (mean $\pm$ standard deviation) of all the methods on 5 simulated datasets. We study the F1 score for the summary causal graph (on the left) and the window-based causal graph (on the right). Results  are computed for various time grid sizes, from 125 to 4000. A log-scale is used in abscissa.}"
              "\n\label{fig:simu:global}\n\end{figure}")
    f_o.close()

