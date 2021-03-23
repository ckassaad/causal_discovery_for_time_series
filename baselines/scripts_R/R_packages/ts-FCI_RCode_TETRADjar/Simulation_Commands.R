# This file contains the commands to run *all* the simulations in "On Causal
# Simdiscovery from Time-Series Data using tsFCI"
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


Simulation_Commands <- function() {
  # simulations in infinite sample size
  # tsfci, fci and Granger

  simInfSamp(0.25, 4, 100)
  simInfSamp(0.5, 4, 100)


  # simulations with tscfci

  # for continuous data 
  simDataCont(0.25, "tscfci", 4, 0.01, 100)
  simDataCont(0.5, "tscfci", 4, 0.01, 100)

  # for discrete data
  simDataDisc(q=0.25, "tscfci", 4, 0.01, 100)
  simDataDisc(q=0.5, "tscfci", 4, 0.01, 100)


  # simulations with cfci

  # for continuous data 
  simDataCont(0.25, "cfci", 4, 0.01, 100)
  simDataCont(0.5, "cfci", 4, 0.01, 100)

  # for discrete data 
  simDataDisc(q=0.25, "cfci", 4, 0.01, 100)
  simDataDisc(q=0.5, "cfci", 4, 0.01, 100)

}
