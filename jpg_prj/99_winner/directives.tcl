############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
############################################################
set_directive_pipeline "PerformIDCT/PerformIDCT_horizontal_loop"
set_directive_array_partition -type complete -factor 8 -dim 2 "DecodeSingleBlock" arrayBlock
set_directive_pipeline "IsInHuffmanCodes"
set_directive_unroll -factor 4 "IsInHuffmanCodes/IsinHuffmanCodes_loop"
