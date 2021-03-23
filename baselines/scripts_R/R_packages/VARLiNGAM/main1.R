main1 <- function(data='macro', path='./', nlags=1, boot_sd = FALSE,
                  boot_irf = FALSE, subsamples = FALSE) {

  # reproduces the results of the paper
  #    'Causal Inference by Independent Component Analysis:
  #                Theory and Applications'
  # by A. Moneta, D. Entner, P.O. Hoyer, and A. Coad

  # INPUT
  # data: 'micro' (firm growth) or 'macro' (monetary polidy)
  # path: path to the data file including file name, f.ex. "myPath/theFile.txt"
  # if data == 'micro' 
  #    nlags = 1 or 2, specifiying the time lags included in the model
  # if data == 'macro'
  #    boot_sd: boolean, if TRUE bootstrap standard deviations of coefficients
  #    boot_irf: boolean, if TRUE bootstrap confidence bands for the irf's
  #    subsamples: boolean, if TRUE do the robustness analysis on subsamples

  sourceDir("./", FALSE)
  sourceDir("lingam/code", FALSE)

  if (data=='macro') {
    main_MonetaryPolicy(path, boot_sd, boot_irf, subsamples)
  }

  if (data=='micro') {
    main_FirmGrowth(path, nlags)
  }

  if (data=='own') {
    print('Loaded files into working space.')
  }

  if (data!='macro' && data!='micro' && data!='own') {
    print("First input must either be 'macro', 'micro' or 'own'.\n")
  }

}
