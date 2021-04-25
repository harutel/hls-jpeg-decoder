############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
############################################################
open_project jpg_prj
set_top JpegDecodeHW
add_files loadjpg.cpp
add_files loadjpg.h
add_files -tb main.cpp
add_files -tb openjpg.cpp
add_files -tb openjpg.h
open_solution "01_no_opt"
set_part {xc7k160tfbg484-1} -tool vivado
create_clock -period 10 -name default
source "./jpg_prj/01_no_opt/directives.tcl"
csim_design -compiler gcc
csynth_design
cosim_design
export_design -format ip_catalog
