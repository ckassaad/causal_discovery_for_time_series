# To start up the R-files
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



# start with 
# > source('start_up.R')
# > start_up()
setwd('/home/kassaad/Documents/Codes/R - codes/ts-FCI/RCode_TETRADjar_tsFCI/RCode_TETRADjar')

start_up <- function() {

  source('dconnected.R')
  source('genData.R')
  source('main_tetrad_fci.R')
  source('plot_timeseries.R')
  source('Plotting_Commands_Barplots.R')
  source('plot_ts_pag.R')
  source('realData_tsfci.R')
  source('scores.R')
  source('Simulation_Commands.R')
  source('Simulations_data_cont.R')
  source('Simulations_data_disc.R')
  source('Simulations_graph.R')
  source('Tetrad_R_interact.R')
  source('ts_functions.R')

}

start_up()

