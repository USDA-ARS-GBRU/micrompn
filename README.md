[![Python package](https://github.com/USDA-ARS-GBRU/micrompn/actions/workflows/python-package.yml/badge.svg)](https://github.com/USDA-ARS-GBRU/micrompn/actions/workflows/python-package.yml)

# MicroMPN: Software for automating most probable number (MPN) estimates from laboratory microplates


## Authors

* Adam Rivers, United States Department of Agriculture, Agricultural Research Service
* Karla Franco Meléndez, United States Department of Agriculture, Agricultural Research Service


## About MPN

Most probable number (MPN) is a  way of estimating how many viable bacteria are in a sample. 
A series of dilutions are made with replication and the  number of  tubes  or wells with bacterial growth is recorded.
This presence/absence and dilution information is used to calculate the most probable number of bacteria in a starting sample.

With liquid handling robots and plate reader this method can be baster than conventional plate counting or 
spotting for some applications. 


With liquid handling robots and plate reader this method can be faster than conventional plate counting or 
spotting on agar for some applications. 



## About MicroMPN

This python package has a command line program processing data from a plate reader. Users provide two types of files: 
A plate layout file specified in Wellmap format and a  CSV plate reader file with columns for the plate, 
the well or row and column and the optical value used for determining growth.

Developers of Python Packages who need an MPN calculation may want to use MicroMPN's mpn estimation function.
This is to my knowledge the only python implementation of MPN estimation in the Python Package Index.
This function is derived from Martine Ferguson and John Ihrie's MPN R package:

Martine Ferguson and John Ihrie. (2019). MPN: Most Probable Number and 
Other Microbial Enumeration Techniques. R package version 0.3.0. <https://CRAN.R-project.org/package=MPN> .



## Installing


The package can be installed from Github 

```bash
  git clone git@github.com:USDA-ARS-GBRU/micrompn.git
  cd micrompn
  pip install .
```


## Getting Started

MicroMPN uses two types of files, a file specifying the layout of the plate and a file containing the data from the plate reader.

1. ["Wellmap"](https://wellmap.readthedocs.io/en/latest/index.html)

- A Wellmap TOML format is a simple way of specifying plate layouts. The Wellmap site has extensive documentation about how to specify the layout. 
   
An example of the format is below. Features can be specifies by rows columns or blocks:

```
[col]
1.dilution = 1e00
2.dilution = 1e-01
3.dilution = 1e-02
4.dilution = 1e-03
5.dilution = 1e-04
6.dilution = 1e-05
7.dilution = 1e-06
8.dilution = 1e-07
9.dilution = 1e-08
10.dilution = 1e-09
11.dilution = 1e-10
12.dilution = 1e-11

[row]
A.replicate = 1
B.replicate = 2
C.replicate = 3
D.replicate = 4
E.replicate = 5
F.replicate = 6
G.replicate = 7
H.replicate = 8
[block.12x8.'A1']
sample = 1
```

Example of the above TOML layout visualized graphically with the Wellmap command "wellmap.show("microplate.toml")

![wellmap replicate](https://github.com/USDA-ARS-GBRU/micrompn/assets/68250738/59f471f0-74e5-44f3-9518-4b287d439745)

![wellmap dilution](https://github.com/USDA-ARS-GBRU/micrompn/assets/68250738/5a1d4203-f26e-46af-8d8c-24580c84e699)


2. Plate reader CSV files

- Plate reader data files vary by instrument but most allow the output of data in tabular format.

- The named columns needed in a CSV are:

    * a plate column
    * a well column OR plate column AND plate row columns
    * an optical value column used for the determination of growth

```
plate_unique, plate_id, plate_well, rfu
plate_0, RFP_1_plate_1_shaking, A1, 27.081
plate_0, RFP_1_plate_1_shaking, A2, 22.001
plate_0, RFP_1_plate_1_shaking, A3, 12.949
plate_0, RFP_1_plate_1_shaking, A4, 10.328
plate_0, RFP_1_plate_1_shaking, A5, 9.264
plate_0, RFP_1_plate_1_shaking, A6, 10.017
plate_0, RFP_1_plate_1_shaking, A7, 9.373
plate_0, RFP_1_plate_1_shaking, A8, 9.049
plate_0, RFP_1_plate_1_shaking, A9, 3.78

```

- The names of these columns are specified in the command-line input "micrompn".

- After running MicroMPN the user will be provided with a CSV containing MPN estimates for each plate and sample

```
plate sample           mpn           mpn_adj         upper         lower
0    plate_0      0  1.005445e+08  9.220112e+07  2.295985e+08  2.295985e+08
1    plate_1      1  1.124383e+08  1.029753e+08  2.609546e+08  2.609546e+08
2    plate_2      2  3.388299e+07  2.797877e+07  6.811641e+07  6.811641e+07
3    plate_3      3  7.636579e+06  6.877494e+06  1.828966e+07  1.828966e+07
4    plate_4      4  7.515884e+05  6.752848e+05  1.195778e+06  1.195778e+06
5    plate_5      5  1.032498e+04  9.468444e+03  1.668272e+04  1.668272e+04
```

- The output contains the MPN, an MPN value corrected for bias due to the number of tubes used and the concentration and an upper and lower bound of the estimate. MicroMPN uses the 95% confidence bound estimation form Jarvis et al. 2010.

Jarvis B, Wilrich C, Wilrich P-T (2010). "Reconsideration of the derivation of Most Probable Numbers, 
their standard deviations, confidence bounds and rarity values." Journal of Applied Microbiology, 1
09, 1660-1667. <https://doi.org/10.1111/j.1365-2672.2010.04792.x>

## Usage

```
MicroMPN: Software to estimate Most Probable Number (MPN) bacterial abundance from
microplates

options:
  -h, --help            show this help message and exit
  --wellmap WELLMAP     A TOML file with plate layout specified in wellmap format
  --data DATA           A csv file or a directory containing csv files with the plate
                        name, optical value, and well or row and column data
  --cutoff CUTOFF       The value from the plate reader above which a well is
                        classified as positive
  --outfile OUTFILE     The file path and name for the results
  --plate_name PLATE_NAME
                        The name of the column containing the plate identifier in the
                        data file
  --value_name VALUE_NAME
                        The name of the column containing the optical signal column
                        in the data file
  --well_name WELL_NAME
                        The name of the column containing the well identifier in the
                        data file
  --col_name COL_NAME   The name of the column containing the plate column identifier
                        in the data file
  --row_name ROW_NAME   The name of the column containing the plate row identifier in
                        the data file
  --zero_padded         If present, the well value in the data file is treated as
                        zero-padded, e.g. A01
  --trim_positives      If present, the list of positive wells will be trimmed to the
                        most dilute all positive dilution and the least dilute all
                        negative dilution. This helps if early dilutions have improbable
                        wells.
  --version, -v         show program's version number and exit
  --logfile LOGFILE, -l LOGFILE
```

The command line options are listed above.

## Example


```bash
micrompn --wellmap micrompn/data/example1_mapfile.toml \
         --data micrompn/data/example1_plate_data.csv \
         --well_name plate_well \
         --plate_name plate_unique \
         --outfile test-output-cutoff6-trimmed.csv \
         --trim_positive \
         --cutoff 6
```


