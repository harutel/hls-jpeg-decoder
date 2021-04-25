# High Level Synthesis of JPEG Decoder

A jpeg decoder implemented in C. The source code has been modified and optimized to be used in Vivado High-Level synthesis tool, which is able to generate RTL designs without the need of writing vhdl or verilog code at all.

## Installation
The operating system which has been used during the development is Ubuntu 16.04 LTS 64-bit and 
This project is implemented in [Vivado Design Suite 2017.2](https://www.xilinx.com/support/download/index.html/content/xilinx/en/downloadNav/vivado-design-tools/archive.html) particularly the HLS Tool.
The only operating system which has been used during development and has been extensively tested is Ubuntu 16.04 LTS 64-bit. Although it's expected to work flawlessly in newer Vivado editions and Ubuntu releases.

## Target Device
Xilinx Kintex 7 xc7k160t fbg484 -1 (xc7k160t:fbg484:-1)

# Usage
To decompress a .jpg file just replace the following two lines in the main function and run the project.
```
char jpg_in[] = "PATH_OF_JPG_FILE_TO_OPEN";
char bmp_out[] = "PATH_OF_BMP_OUTPUT_FILE";
```
Please read Vivado Tutorials for more information how to use the tool.

# Limitations
It's very important to be considered  is the fact that

## Contributing
Since this project is mainly developed as a part of a MSc Thesis the quality of the code is very poor. Therefore it's highly recommended to be used only for experimental use.
Any contribution to this projects is welcome and appreciated.

## Cite
[High Level Synthesis of JPEG Decoder on FPGA](http://hdl.handle.net/10889/10581)
Source code is based on [Jpeg (.jpg) Image Source Decoder](https://xbdev.net/image_formats/jpeg/jpeg_decoder_source) demo by bkenwright@xbdev.net

## License
[MIT](https://choosealicense.com/licenses/mit/)
