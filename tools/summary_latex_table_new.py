if __name__ == "__main__":
    measure_name = "F-Score"  # Precision Recall F-Score

    benchmark = 'fmri'

    name_of_method = {'GrangerPW': 'PWGC',
                        'GrangerMV': 'MVGC',
                        'TCDF': 'TCDF',
                        'PCMCICMIknn': 'PCMCImi',
                        'PCMCIParCorr': 'PCMCIpc',
                        'PCTMI': 'PCTMI',
                        'oCSE': 'oCSE',
                        'tsFCI': 'tsFCI',
                        'FCITMI': 'FCITMI',
                        'VarLiNGAM': 'VarLiNGAM',
                        'TiMINO': 'TiMINO'}

    path_output = "../experiments/performance_average/latex_format/tables" + "/" + measure_name + ".txt"
    f_o = open(path_output, "w")
    f_o.write("\\begin{table*}[ht!]\n\t\centering\n")
    f_o.write("\t\\begin{tabular}{c|c|c|c|c|c}\n & \\multicolumn{1}{c|}{Summary} &  \\multicolumn{3}{c}{Summary Details} \\\\ \n && Ext. Adj. & Ext. Caus. & Self Caus.  \\\\ \\hline\n")
    structure = benchmark
    measure = ""
    all_methods_list = ['GrangerPW', 'GrangerMV', 'TCDF', 'PCMCICMIknn', 'PCMCIParCorr',
                        'oCSE', 'tsFCI', 'VarLiNGAM', 'TiMINO']
    for method in all_methods_list:
        f_o.write("\t\t"+name_of_method[method])
        for type_of_causes in ["summary_other_and_self", "summary_other", "summary_other", "summary_self"]:
            path_input = "../experiments/performance_average/" + type_of_causes + "_performance_average"
            if type_of_causes == "summary_other_and_self":
                measure_type = "Oriented"
                measure = measure_name + " " + measure_type
            elif type_of_causes == "summary_other":
                measure_type1 = "Adjacent"
                measure_type2 = "Oriented"
                if measure == measure_name + " " + measure_type1:
                    measure = measure_name + " " + measure_type2
                else:
                    measure = measure_name + " " + measure_type1
            elif type_of_causes == "summary_self":
                measure = measure_name
            else:
                all_methods_list = None
                methods_treats_hidden = None

            files_input_name = method + "_" + structure
            print(files_input_name)
            try:
                f_i = open(str(path_input) + "/" + str(files_input_name), "r")
                for line_i in f_i:
                    if measure in line_i:
                        nextLine = next(f_i)
                        print(method, structure, type_of_causes)
                        print(nextLine)
                        nextLine = nextLine.replace("\n", "")
                        nextline_processed = nextLine.split(" ")
                        nextLine = nextLine.replace("+-", "\\pm")
                        measure_mean = "{:.2f}".format(round(float(nextline_processed[0]), 2))
                        measure_std = "{:.2f}".format(round(float(nextline_processed[2]), 2))
                        result_line = "&$"+measure_mean+" \\pm "+measure_std+"$"
                        print(result_line)
                        f_o.write(result_line)
                        break
                f_i.close()
            except FileNotFoundError:
                print(files_input_name)
        f_o.write("\\\\ \n")
    f_o.write("\t\\end{tabular} \\caption{Results for FMRI datasets. We provide the F1 score (mean $\\pm$ standard deviation over the 27 datasets considered)} \n \\end{table*}")
    f_o.close()

