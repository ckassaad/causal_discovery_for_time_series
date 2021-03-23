# This file contains the commands to make the bar-plots of the main results
# in "On Causal discovery from Time-Series Data using tsFCI"
# 
# Copyright (C) 2010 Doris Entner
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see http://www.gnu.org/licenses/


plot_commands <- function() {

  # plots all three scores for all methods and sample sizes and q's
  # figures will be stored to simulation_results/figures/
  
  path_loadfrom <- "./simulation_results/scores/"
#  path_loadfrom <- "./simulation_results"

  # ---------------------------------------------------------------------------
  # discrete data

  nfiles <- 6
  filename <- array("",dim=c(1,nfiles))

  # q = 0.25
  filename[1] <- paste(path_loadfrom, "RESdisc25tscfci4_1.txt", sep="") # tscfci
  filename[2] <- paste(path_loadfrom, "RESdisc25cfci4_1.txt", sep="") # cfci
  filename[3] <- paste(path_loadfrom, "RESinf25_4.txt", sep="") # infinite sample
  # q = 0.5
  filename[4] <- paste(path_loadfrom, "RESdisc50tscfci4_1.txt", sep="") # tscfci
  filename[5] <- paste(path_loadfrom, "RESdisc50cfci4_1.txt", sep="") # cfci
  filename[6] <- paste(path_loadfrom, "RESinf50_4.txt", sep="") # infinite sample

  # load scores
  s <- vector(length=nfiles,mode='list')
  for (i in 1:nfiles) {
    s[[i]] <- readListFromFile(filename[i])
  }

  # summarize scores
  s_sum <- matrix(0,9*nfiles,2)
  cnt <- 0
  for (i in 1:nfiles) {
    for (j in 1:length(s[[i]])) {
      cnt <- cnt+1
      s_sum[cnt,] <- colSums(s[[i]][[j]])/nrow(s[[i]][[j]])
    }
  }

  # reorder: discrete - tsCFCI, FCI - 100, 1000, 10000, inf - q=0.25, q=0.5
  #          ancestor
  #          pairwise
  # take out Granger in infinite sample size!
  # see my papers p.88
  s_DC_d <- s_sum[c(1,10,2,11,3,12,19,21,28,37,29,38,30,39,46,48),]
  s_AN_d <- s_sum[c(4,13,5,14,6,14,22,24,31,40,32,41,33,42,49,51),]
  s_PW_d <- s_sum[c(7,16,8,17,9,18,25,27,34,43,35,44,36,45,52,54),]

  # convert scores in percentages of right and wrong decisions
  sp_DC_d <- matrix(0,16,2)
  sp_AN_d <- matrix(0,16,2)
  sp_PW_d <- matrix(0,16,2)
  for (i in 1:16) {
    sp_DC_d[i,1] <- s_DC_d[i,1]*s_DC_d[i,2] # percentage of right decisions
    sp_DC_d[i,2] <- s_DC_d[i,1]-sp_DC_d[i,1] # percentage of wrong decisions
    sp_AN_d[i,1] <- s_AN_d[i,1]*s_AN_d[i,2]
    sp_AN_d[i,2] <- s_AN_d[i,1]-sp_AN_d[i,1]
    sp_PW_d[i,1] <- s_PW_d[i,1]*s_PW_d[i,2]
    sp_PW_d[i,2] <- s_PW_d[i,1]-sp_PW_d[i,1]
  }

  sp_DC_d <- sp_DC_d*100
  sp_AN_d <- sp_AN_d*100
  sp_PW_d <- sp_PW_d*100

  # ---------------------------------------------------------------------------
  # continuous data
  nfiles <- 6
  filename <- array("",dim=c(1,nfiles))

  # q = 0.25
  filename[1] <- paste(path_loadfrom, "REScont25tscfci4_1.txt", sep="") # tscfci
  filename[2] <- paste(path_loadfrom, "REScont25cfci4_1.txt", sep="") # cfci
  filename[3] <- paste(path_loadfrom, "RESinf25_4.txt", sep="") # infinite sample
  # q = 0.5
  filename[4] <- paste(path_loadfrom, "REScont50tscfci4_1.txt", sep="") # tscfci
  filename[5] <- paste(path_loadfrom, "REScont50cfci4_1.txt", sep="") # cfci
  filename[6] <- paste(path_loadfrom, "RESinf50_4.txt", sep="") # infinite sample

  # load scores
  s <- vector(length=nfiles,mode='list')
  for (i in 1:nfiles) {
    s[[i]] <- readListFromFile(filename[i])
  }

  # summarize scores
  s_sum <- matrix(0,9*nfiles,2)
  cnt <- 0
  for (i in 1:nfiles) {
    for (j in 1:length(s[[i]])) {
      cnt <- cnt+1
      s_sum[cnt,] <- colSums(s[[i]][[j]])/nrow(s[[i]][[j]])
    }
  }

  # add Granger, PSI and Group Lasso scores (from matlab)
  s2 <- matrix(0,30,2)
  # --- q = 0.25
  # score_GraDC =
  s2[1,] <- c(1.0000, 0.9337)
  s2[2,] <- c(1.0000, 0.9433)
  s2[3,] <- c(1.0000, 0.9216)
  # score_GraAN =
  s2[4,] <- c(1.0000, 0.8928)
  s2[5,] <- c(1.0000, 0.9196)
  s2[6,] <- c(1.0000, 0.9002)
  # score_GraPW =
  s2[7,] <- c(1.0000, 0.9039)
  s2[8,] <- c(1.0000, 0.9535)
  s2[9,] <- c(1.0000, 0.9416)
  # score_PsiPW =
  s2[10,] <- c(0.8431, 0.7671)
  s2[11,] <- c(0.7018, 0.8096)
  s2[12,] <- c(0.6493, 0.8091)
  # score_GLaPW =
  s2[13,] <- c(1.0000, 0.8648)
  s2[14,] <- c(1.0000, 0.9031)
  s2[15,] <- c(1.0000, 0.9097)

  # --- q = 0.5
  # score_GraDC50 =
  s2[16,] <- c(1.0000, 0.7874)
  s2[17,] <- c(1.0000, 0.8418)
  s2[18,] <- c(1.0000, 0.8331)
  # score_GraAN50 =
  s2[19,] <- c(1.0000, 0.8498)
  s2[20,] <- c(1.0000, 0.9151)
  s2[21,] <- c(1.0000, 0.8964)
  # score_GraPW50 =
  s2[22,] <- c(1.0000, 0.8470)
  s2[23,] <- c(1.0000, 0.9465)
  s2[24,] <- c(1.0000, 0.9383)
  # score_PsiPW50 =
  s2[25,] <- c(0.7692, 0.5408)
  s2[26,] <- c(0.5980, 0.7646)
  s2[27,] <- c(0.5419, 0.8439)
  # score_GLaPW50 =
  s2[28,] <- c(1.0000, 0.8234)
  s2[29,] <- c(1.0000, 0.8386)
  s2[30,] <- c(1.0000, 0.8489)

  s_sum <- rbind(s_sum,s2)


  # reorder: discrete - tsCFCI, FCI, Granger - 100, 1000, 10000, inf - q=0.25, q=0.5
  #          ancestor - tsCFCI, FCI, Granger
  #          pairwise - tsCFCI, FCI, Granger, PSI, Group Lasso
  # see my papers p.88f
  s_DC_c <- s_sum[c(1,10,55,2,11,56,3,12,57,19,21,20, 28,37,70,29,38,71,30,39,72,46,48,47),]
  s_AN_c <- s_sum[c(4,31,58,5,32,59,6,33,60,22,24,23, 31,40,73,32,41,74,33,42,75,49,51,50),]
  s_PW_c <- s_sum[c( 7,16,61,64,67, 8,17,62,65,68, 9,18,63,66,69,25,27,26, 
                  34,43,76,79,82,35,44,77,80,83,36,45,78,81,84,52,54,53),]

  # convert scores in percentages of right and wrong decisions
  sp_DC_c <- matrix(0,24,2)
  sp_AN_c <- matrix(0,24,2)
  sp_PW_c <- matrix(0,36,2)
  for (i in 1:24) {
    sp_DC_c[i,1] <- s_DC_c[i,1]*s_DC_c[i,2] # percentage of right decisions
    sp_DC_c[i,2] <- s_DC_c[i,1]-sp_DC_c[i,1] # percentage of wrong decisions
    sp_AN_c[i,1] <- s_AN_c[i,1]*s_AN_c[i,2]
    sp_AN_c[i,2] <- s_AN_c[i,1]-sp_AN_c[i,1]
    sp_PW_c[i,1] <- s_PW_c[i,1]*s_PW_c[i,2]
    sp_PW_c[i,2] <- s_PW_c[i,1]-sp_PW_c[i,1]
  }
  for (i in 25:36) {
    sp_PW_c[i,1] <- s_PW_c[i,1]*s_PW_c[i,2]
    sp_PW_c[i,2] <- s_PW_c[i,1]-sp_PW_c[i,1]
  }

  sp_DC_c <- sp_DC_c*100
  sp_AN_c <- sp_AN_c*100
  sp_PW_c <- sp_PW_c*100


  # ------------------------------ plotting -----------------------------------

  # for discrete (DC, AN, PW)
  space_bp1 <- rep(c(0,0.1,0.5,0.1,0.5,0.1,0.5,0.1),2)
  space_bp1[9] <- 1 # length=16

  # for continuous (DC, AN)
  space_bp2 <- rep(c(0,0.1,0.1,0.5,0.1,0.1,0.5,0.1,0.1,0.5,0.1,0.1),2)
  space_bp2[13] <- 1 # length=24

  # for continuous (PW)
  space_bp3 <- rep(c(0,0.1,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.1,0.5,0.1,0.1),2)
  space_bp3[19] <- 1 # length=36


  # ----- DC and AN scores
  sp_DC <- rbind(sp_DC_d,sp_DC_c) # all direct-cause scores
  sp_AN <- rbind(sp_AN_d,sp_AN_c) # all ancestor scores

  space_bp <- c(space_bp1,space_bp2)
  space_bp[17] <- 2

  # font style in plot
  myfont = "Times"; #"Helvetica"

  # --- DC 
  pdf(file="./simulation_results/figures/DC_scores.pdf", width=16, height=6, 
      paper="special",family=myfont)

#  x11(width=16, height=6)
  bp <- plot_mybars_new(sp_DC,space_bp)
  mtext("Direct-cause scores",side=3, at=(length(space_bp)+sum(space_bp)-1)/2, line=1, cex=2)

  add_text_DC_AN(bp)
  dev.off()

  # --- AN
  pdf(file="./simulation_results/figures/AN_scores.pdf", width=16, height=6, 
      paper="special",family=myfont)
  bp <- plot_mybars_new(sp_AN,space_bp)
  mtext("Ancestor scores",side=3, at=(length(space_bp)+sum(space_bp)-1)/2, line=1, cex=2)

  add_text_DC_AN(bp)
  dev.off()

  # ----- PW scores
  sp_PW <- rbind(sp_PW_d,sp_PW_c) # all pairwise scores

  space_bp <- c(space_bp1,space_bp3)
  space_bp[17] <- 2

  pdf(file="./simulation_results/figures/PW_scores.pdf", width=20, height=6, 
      paper="special",family=myfont)
  bp <- plot_mybars_new(sp_PW,space_bp)
  mtext("Pairwise scores",side=3, at=(length(space_bp)+sum(space_bp)-1)/2, line=1, cex=2)

  add_text_PW(bp)
  dev.off()

}


plot_mybars_new <- function(scores, spaces) {

  par(mar=c(10.5,3.6,3.6,1.1))

  mygreen <- "#ccffcc"

  xlimits <- c(1,length(spaces)+sum(spaces)-1)
  ylimits <- c(-0.1,100.5)

  bp_pos <- barplot(t(scores), xlim=xlimits, ylim=ylimits,
            space=spaces, density=c(1,1), angle=c(135,45),
            col=c(mygreen,"red"),yaxt="n") 
  axis(2, at=c(0,50,100),labels=c(0,50,100),cex.axis=2, las=2)

  grid(nx = 0, ny = 10, col = "gray", lty = "solid", lwd=0.2)

  # to get grid behind the bars - just plot everything again...
  oldpar <- par(bg='white')
  par(new=T) 
  barplot(t(scores), xlim=xlimits, ylim=ylimits, 
            space=spaces,col=c(mygreen,"white"),yaxt="n") 
  par(new=T)
  bp_pos <- barplot(t(scores), xlim=xlimits, ylim=ylimits,
            space=spaces, density=c(1,16), angle=c(135,45),
            col=c(mygreen,"red"),yaxt="n") 
  par(oldpar)

  abline(h=c(0,100))

  bp_pos

}


add_text_DC_AN <- function(bp) {

  # add text to plot
  # in bp are the positions of the bars

  sample <- rep(c('T=100','T=1000','T=10000',expression(T %->% infinity)),4)
  pos_sample <- c()
  for (i in 1:8*2-1) {
    pos_sample <- c(pos_sample, (bp[i]+bp[i+1])/2)
  }
  for (i in seq(18,40,by=3)) {
    pos_sample <- c(pos_sample, bp[i])
  }

  sparsity <- rep(c("q=0.25","q=0.5"),2)
  pos_sparsity <- c((bp[4]+bp[5])/2, (bp[12]+bp[13])/2, (bp[22]+bp[23])/2, (bp[34]+bp[35])/2)

  datatype <- c("binary data", "continuous-valued data")
  pos_datatype <- c((bp[8]+bp[9])/2, (bp[28]+bp[29])/2)

  # side: on which side of the plot (1=bottom, 2=left, 3=top, 4=right)
  mtext(sample, side=1, at=pos_sample, line=1, cex=1)
  mtext(sparsity, side=1, at=pos_sparsity, line=3, cex=1.5)
  mtext(datatype, side=1, at=pos_datatype, line=5, cex=2)

  algnames <- rep(c("tsFCI", "FCI"),8)
  algnames <- c(algnames, rep(c("tsFCI", "FCI", "Granger"),8))

  text(bp, 5, algnames, srt=90, pos=4, cex=1.3, offset = 0.08, font=2) #cex=1)

}


add_text_PW <- function(bp) {

  # add text to plot
  # in bp are the positions of the bars
  
  sample <- rep(c('T=100','T=1000','T=10000',expression(T %->% infinity)),4)
  pos_sample <- c()
  for (i in 1:8*2-1) {
    pos_sample <- c(pos_sample, (bp[i]+bp[i+1])/2)
  }
  for (i in c(19,24,29,33,37,42,47,51)) {
    pos_sample <- c(pos_sample, bp[i])
  }

  sparsity <- rep(c("q=0.25","q=0.5"),2)
  pos_sparsity <- c((bp[4]+bp[5])/2, (bp[12]+bp[13])/2, bp[26], bp[44])

  datatype <- c("binary data", "continuous-valued data")
  pos_datatype <- c((bp[8]+bp[9])/2, (bp[34]+bp[35])/2)

  # side: on which side of the plot (1=bottom, 2=left, 3=top, 4=right)
  mtext(sample, side=1, at=pos_sample, line=1, cex=1)
  mtext(sparsity, side=1, at=pos_sparsity, line=3, cex=1.2)
  mtext(datatype, side=1, at=pos_datatype, line=5, cex=2)

  algnames <- rep(c("tsFCI", "FCI"),8)
  algnames1 <- rep(c("tsFCI", "FCI", "Granger", "PSI", "Group Lasso"),3)
  algnames1 <- c(algnames1,c("tsFCI", "FCI", "Granger"))
  algnames1 <- rep(algnames1,2)
  algnames <- c(algnames,algnames1)

  text(bp, 5, algnames, srt=90, pos=4, cex=1.3, offset = 0.08, font=2) #cex=1)

}
