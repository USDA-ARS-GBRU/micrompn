"""MicroMPN: a Python command line program for automating microbiology most probabable 
   number (MPN) estimates in microtiter paltes

"""
import os
import argparse
import logging
from typing import Tuple, List, Dict

import numpy as np
import pandas as pd
import wellmap
import micrompn.mpn


def myparser() -> argparse.ArgumentParser:
    """Parse the input arguments for micrompn 

    :return: _description_
    :rtype: argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser(
    description='micrompn: Software to estimate Most Probable Number (MPN) bacterial abundance from microtiter plates')
    parser.add_argument('--wellmap', type=str, required=True, help='A TOML file  with plate layout speficied in wellmap format ')
    parser.add_argument('--data', type=str, required=True, help='A csv file or a directory contiining csv files with the plate name, optical value, and well or row and column data')
    parser.add_argument('--cutoff', type=float, required=True, help='The value from the plate reader above which a well is classfied as positive')
    parser.add_argument('--outfile', type=str, required=True, help='The file path and name for the results')
    parser.add_argument('--plate_name', type=str, default='plate', help='The name of the column containing the plate identifier in the data file')
    parser.add_argument('--value_name', type=str, default='rfu', help='The name of the column containing the optical signal column in the data file')
    parser.add_argument('--well_name', type=str, required=False, help='The name of the column containing the well identifier in the data file')
    parser.add_argument('--col_name', type=str, required=False, help='The name of the column containing the plate column identifier in the data file')
    parser.add_argument('--row_name', type=str, required=False, help='The name of the column containing the plate row identifier in the data file')   
    parser.add_argument('--zero_padded', action='store_true', help='if present the well value in the data file is treated as zero-padded, e.g. A01')
    parser.add_argument('--version', '-v',  action='version', version="%(prog)s (" + micrompn.__version__ + ")")
    parser.add_argument('--logfile', '-l', type=str, default= "micrompn.log")
    return parser
    
def _logger_setup(logfile: str) -> logging.Logger:
    """Set up logger for micrompn

    :param logfile: The path for the logger
    :type logfile: str
    :raises e: A generic error if logger could not be set up
    :return: Returns a confgured logger
    :rtype: logging.logger
    """
    try:
        # create logger
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        # create console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create file handler

        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
        # tell the handler to use this format
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(fh)
        logger.addHandler(ch)
        return logger
    except Exception as e:
        print("An error occurred setting up logging")
        raise e


def split_dataframe(df: pd.DataFrame) -> Tuple[List[pd.DataFrame], pd.Series]:
    """ The function preaks a dataframe of plate data and layout data and into a Tuple cotnining a list of 
        dataframes split by plates and samples and a pd.Series of unique samples

    :param df: dataframe containing a sample column
    :type df: pd.DataFrame
    :return:  list of dataframse split by sample and a Pandas Series with the samples names 
    :rtype: Tuple[List[pd.DataFrame], pd.Series]
    """
    df['plate_sample'] = df["plate"] + '_' + df['sample'].astype(str)
    unique_samples = df['plate_sample'].unique() # Get unique values in the 'plate_sample' column
    list_of_dfs = []
    for sample in unique_samples:
        sample_df = df[df['plate_sample'] == sample] # Filter the dataframe by the unique sample
        list_of_dfs.append(sample_df) # Add the filtered dataframe to the list
    return list_of_dfs, unique_samples

def count_above_cutoff(df: pd.DataFrame, cutoff: float) -> np.array:
    """ Counts how many wells are above a cutoff value at each dilution

    :param df: The input dataframe containing dilution values
    :type df: pd.DataFrame
    :param cutoff: The value above which a well is considered positive
    :type cutoff: float
    :return: An array of the counts of positive wells at each dilution
    :rtype: np.array
    """
    grouped = df.groupby('dilution')
    counts = []
    for name, group in grouped:
        above_cutoff = np.sum(group['value'] > cutoff)
        counts.append(above_cutoff)
    counts.reverse()
    return np.array(counts)


def import_plate(platefile: str,
                 plate_name: str,
                 well_name: str,
                 col_name: str,
                 row_name: str,
                 value_name: str,
                 zero_padded: bool) -> pd.DataFrame:
    """ import a csv of column data from a plate reader and convert it to a uniform 
        format for merging with the  plate layout file. This works for 96 and 384 well plates but not 1536 wells

    :param platefile: A string of the path to the file.
    :type platefile: str
    :param plate_name: A column label with a unique plate name.
    :type plate_name: str
    :param well_name: A column label for the well name (e.g. A1 or A01), optional if col and row are provided.
    :type well_name: str
    :param col_name: A column label for the column name (e.g. 12) rownames should br proveded too. optional if wellname is provided.
    :type col_name: str
    :param row_name: A column label for the row name (e.g. D) colnames should br proveded too. optional if wellname is provided.
    :type row_name: str
    :param value_name: A column label for the optical measurment used to determine Growth in the well
    :type value_name: str
    :param zero_padded: Does well have leading zero in the column number, e.g.  A01 instead of A1
    :type zero_padded: bool
    :return: A uniform dataframe for merging with the plate layout dataframe.
    :rtype: pd.DataFrame
    """
    assert(well_name or (col_name and row_name)),"you need to specify the columns to use for microplate well or microbplate column and row."
    try:
        plate = pd.read_csv(platefile)
        dtype_dict = {'plate' : None, 'well' : None, 'col' : None, 'row' : None, 'value' : None }
        cleanplate = pd.DataFrame(columns=dtype_dict.keys()).astype(dtype_dict)
        cleanplate['plate'] = plate[plate_name]
        cleanplate['value'] = plate[value_name]
        if col_name and row_name:
            cleanplate['col'] = plate[col_name]
            cleanplate['row'] = plate[row_name]
            cleanplate['well'] =  plate[row_name] +  plate[col_name].str()
        elif well_name:
            if zero_padded:
                cleanplate['well'] = plate[well_name].astype(str).str[0] + plate[well_name].astype(str).str[1:].int().str()
            else:
                cleanplate['well'] = plate[well_name]    
            cleanplate['row'] = plate[well_name].str[0]
            cleanplate['col'] = plate[well_name].str.slice(start=1).str.zfill(2)
        return cleanplate
    except:
        raise Exception
        print("could not parse input plate data")
    
def process_plate(wellmap_file: str, plate_df: pd.DataFrame, args: argparse.Namespace ) -> pd.DataFrame:
    """Read plate data and return a dataframe of MPN estimates and statistics for each sample

    :param wellmap_file: Path to the Wellmap layout TOML file
    :type wellmap_file: str
    :param plate_df: a pandas DataFrame of plae reader data from one or more plates
    :type plate_df: pd.DataFrame
    :param args: The input parameter arguments
    :type args: argparse.Namespace
    :return: A data frame of plate, sample, MPN and MPN statistics
    :rtype: pd.DataFrame
    """
    dtype_dict = {'plate': str, 'sample': str, 'mpn': float, 'mpn_adj': float, 'upper': float, 'lower': float}
    temp_df= pd.DataFrame(columns=dtype_dict.keys()).astype(dtype_dict)
    layout = wellmap.load(toml_path=wellmap_file)
    data = layout.merge(right=plate_df, on = "well", how='left')
    df_list, unique_samples = split_dataframe(data)
    for df, sample_no in zip(df_list, unique_samples):
        positive = count_above_cutoff(df, args.cutoff)
        dilution_list = (df['dilution'].unique())
        tubes = np.repeat(a = len(df['replicate'].unique()), repeats=len(dilution_list))
        mobj = micrompn.mpn(positive, tubes, dilution_list) 
        temp_df.loc[len(temp_df)] = [list(df['plate'])[0], sample_no.split('_')[1], mobj['MPN'], mobj['MPN_adj'], mobj['UB'], mobj['UB']]
    return temp_df


def is_file_or_directory(path: str) -> bool:
    """check if a string is a path or a directory

    :param path: the file or director
    :type path: str
    :raises FileNotFoundError: Error is not a file or directoty
    :return: true if a file, false if a directory
    :rtype: bool
    """
    if os.path.isfile(path):
        return True
    elif os.path.isdir(path):
        return False
    else:
        raise FileNotFoundError(f"{path} not found")

def main(arglist: list = None) -> None:
    """An entry point to run the full micrompn workflow

    :param arglist: A list of arguments, used for testing, defaults to None
    :type arglist: list, optional
    """
    # Set up logging
    parser = myparser()
    args = parser.parse_args(arglist)
    logger = _logger_setup(args.logfile)
    # set up pandas datafreame 
    dtype_dict = {'plate': str, 'sample': str, 'mpn': float, 'mpn_adj': float, 'upper': float, 'lower': float}
    temp_df= pd.DataFrame(columns=dtype_dict.keys()).astype(dtype_dict)
    if is_file_or_directory(args.data):
        clean_plate = import_plate(platefile=args.data,
                                   plate_name=args.plate_name,
                                   well_name=args.well_name,
                                   col_name=args.col_name,
                                   row_name=args.row_name,
                                   value_name=args.value_name,
                                   zero_padded=args.zero_padded)
        newplate_df = process_plate(args.wellmap, clean_plate, args=args)
        temp_df = pd.concat([temp_df, newplate_df], ignore_index = True)
    else:
        for plate_file in os.listdir(args.datadir):
            clean_plate = import_plate(platefile=plate_file,
                                       plate_name=args.plate_name,
                                       well_name=args.well_name,
                                       col_name=args.col_name,
                                       row_name=args.row_name,
                                       value_name=args.value_name,
                                       zero_padded=args.zero_padded)
            newplate_df = process_plate(args.wellmap, clean_plate, args=args)
            temp_df = pd.concat([temp_df, newplate_df], ignore_index = True)
    print(temp_df)
    temp_df.to_csv(args.outfile)