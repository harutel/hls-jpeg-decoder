############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
############################################################
open_project jpg_prj
set_top JpegDecodeHW
add_files src/loadjpg.cpp
add_files -tb src/main.cpp
add_files -tb src/openjpg.cpp
open_solution "99_winner"
set_part {xc7k160tfbg484-1}
create_clock -period 10 -name default
source "./jpg_prj/99_winner/directives.tcl"
csim_design -compiler gcc
csynth_design
cosim_design
export_design -format ip_catalog
