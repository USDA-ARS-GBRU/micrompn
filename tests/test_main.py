import argparse

import pytest
import numpy as np
import pandas as pd


from micrompn import (
    myparser,
    split_dataframe,
    count_above_cutoff,
    trim_counts,
    import_plate,
    process_plate,
    is_file_or_directory,
    main
)

def test_myparser():
    parser = myparser()

    assert isinstance(parser, argparse.ArgumentParser)
    assert parser.description == 'MicroMPN: Software to estimate Most Probable Number (MPN) bacterial abundance from microplates'
    assert any(action.option_strings[0] == '--wellmap' for action in parser._actions)
    assert any(action.option_strings[0] == '--data' for action in parser._actions)
    assert any(action.option_strings[0] == '--cutoff' for action in parser._actions)
    assert any(action.option_strings[0] == '--outfile' for action in parser._actions)
    assert any(action.option_strings[0] == '--plate_name' for action in parser._actions)
    assert any(action.option_strings[0] == '--value_name' for action in parser._actions)
    assert any(action.option_strings[0] == '--well_name' for action in parser._actions)
    assert any(action.option_strings[0] == '--col_name' for action in parser._actions)
    assert any(action.option_strings[0] == '--row_name' for action in parser._actions)
    assert any(action.option_strings[0] == '--zero_padded' for action in parser._actions)
    assert any(action.option_strings[0] == '--trim_positives' for action in parser._actions)
    assert any(action.option_strings[0] == '--version' for action in parser._actions)
    assert any(action.option_strings[0] == '--logfile' for action in parser._actions)
    assert parser._option_string_actions['--version'].dest == 'version'
    assert parser._option_string_actions['--logfile'].dest == 'logfile'

def test_split_dataframe():
    # Create a sample dataframe


    df = pd.DataFrame({
        'plate': ['Plate1', 'Plate1', 'Plate2', 'Plate2'],
        'sample': ['Sample1', 'Sample1', 'Sample2', 'Sample2'],
        'value': [1., 2., 3., 4.]
    })

    # Call the function to split the dataframe
    result_dfs, result_unique_samples = split_dataframe(df)

    # Verify the type and length of the returned values
    assert isinstance(result_dfs, list)
    assert isinstance(result_unique_samples, np.ndarray)
    assert len(result_dfs) == len(result_unique_samples)

    # Verify the content of the split dataframes
    expected_dfs = [
        df[df['plate_sample'] == 'Plate1_Sample1'],
        df[df['plate_sample'] == 'Plate2_Sample2'],
    ]
    for expected_df, result_df in zip(expected_dfs, result_dfs):
        assert result_df.equals(expected_df)

    # Verify the content of the unique samples series
    expected_unique_samples = pd.Series(['Plate1_Sample1', 'Plate2_Sample2'])
    assert pd.Series(result_unique_samples).equals(expected_unique_samples)

def test_count_above_cutoff():
    # Create a sample dataframe
    df = pd.DataFrame({
        'dilution': [1.0, 1.0, 1.0, 0.1, 0.1, 0.1],
        'value': [5.5, 4.5, 3.5, 3.5, 2.5, 1.5]
    })
    cutoff = 3.0

    # Call the function to count positive wells above the cutoff
    result = count_above_cutoff(df, cutoff)

    # Verify the type and length of the returned array
    assert isinstance(result, np.ndarray)
    assert len(result) == len(df['dilution'].unique())

    # Verify the content of the counts array
    expected_counts = np.array([3,1])
    assert np.array_equal(result, expected_counts)

def test_trim_counts():
    # Test case 1: Trim all positive and negative dilutions exist
    positive = np.array([10, 10, 5, 1, 0, 0])
    tubes = np.array([10, 10, 10, 10, 10, 10])
    dilution = np.array([1., 0.1, 0.01, 0.001, 0.0001, 0.00001])
    expected_result = (
        np.array([10, 5, 1, 0]),
        np.array([10, 10, 10, 10]),
        np.array([0.1, 0.01, 0.001, 0.0001])
    )

    result = trim_counts(positive, tubes, dilution)
    print(result)
    assert isinstance(result, tuple)
    assert len(result) == 3
    assert np.array_equal(result[0], expected_result[0])
    assert np.array_equal(result[1], expected_result[1])
    assert np.array_equal(result[2], expected_result[2])

    # Test case 2: No all positive or all negative dilutions
    positive = np.array([1, 0, 1, 0])
    tubes = np.array([10, 10, 10, 10])
    dilution = np.array([1, 2, 3, 4])

    result = trim_counts(positive, tubes, dilution)
    assert isinstance(result, tuple)
    assert len(result) == 3
    assert np.array_equal(result[0], positive)
    assert np.array_equal(result[1], tubes)
    assert np.array_equal(result[2], dilution)

    # Test case 3: Empty arrays
    positive = np.array([])
    tubes = np.array([])
    dilution = np.array([])

    result = trim_counts(positive, tubes, dilution)
    assert isinstance(result, tuple)
    assert len(result) == 3
    assert np.array_equal(result[0], positive)
    assert np.array_equal(result[1], tubes)
    assert np.array_equal(result[2], dilution)
