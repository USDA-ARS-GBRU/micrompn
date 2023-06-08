import numpy as np
import pytest
#from scipy.optimize import root_scalar
#from scipy.stats import norm

from micrompn import (
    ptEst_MPN,
    ptEstAdj_MPN,
    logL_MPN,
    logLR_MPN,
    logLRroot_MPN,
    jarvisCI_MPN,
    likeRatioCI_MPN,
    rarity_MPN,
    isMissing,
    checkInputs_mpn,
    mpn,
)

def test_ptEst_MPN():
    # Test case 1: All negative tubes
    positive = np.array([0, 0, 0])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    assert ptEst_MPN(positive, tubes, amount) == 0

    # Test case 2: All positive tubes
    positive = np.array([10, 10, 10])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    assert ptEst_MPN(positive, tubes, amount) == float('inf')

    # Test case 3: General case
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    expected_result = 74.22876396932874
    assert expected_result == pytest.approx(ptEst_MPN(positive, tubes, amount))

def test_ptEstAdj_MPN():
    MPN = 74.22876396932874
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    expected_result = 67.92847966627644

    result = ptEstAdj_MPN(MPN, positive, tubes, amount)

    assert isinstance(result, float)
    assert np.isfinite(result)
    assert np.isclose(result, expected_result)

def test_logL_MPN():
    lambda_ = 74.22876396932874
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    expected_result = -10.254505611153164

    result = logL_MPN(lambda_, positive, tubes, amount)

    assert isinstance(result, float)
    assert np.isclose(result, expected_result)

def test_logLR_MPN():
    lambda_ = 74.22876396932874
    lambda_hat = 73.0
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    expected_result = -0.001695304319877522

    result = logLR_MPN(lambda_, lambda_hat, positive, tubes, amount)

    assert isinstance(result, float)
    assert np.isclose(result, expected_result)

def test_logLRroot_MPN():
    lambda_ = 74.22876396932874
    lambda_hat = 73.0
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    crit_val = 0.95
    expected_result = -0.9516953043198775

    result = logLRroot_MPN(lambda_, lambda_hat, positive, tubes, amount, crit_val)

    assert isinstance(result, float)
    assert np.isclose(result, expected_result)

def test_jarvisCI_MPN():
    MPN = 74.22876396932874
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    conf_level = 0.95

    expected_output = {
        'variance': 902.4173735493894,
        'var_logMPN': 0.1637808007299598,
        'LB': 33.58089913036114,
        'UB': 164.07867398144504
    }

    result = jarvisCI_MPN(MPN, positive, tubes, amount, conf_level)

    assert isinstance(result, dict)
    assert result == expected_output

def test_likeRatioCI_MPN():
    MPN = 74.22876396932874
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])
    conf_level = 0.95

    expected_output = {
        'LB': 32.43393676359959,
        'UB': 152.3661892483937
    }

    result = likeRatioCI_MPN(MPN, positive, tubes, amount, conf_level)

    assert isinstance(result, dict)
    assert np.isclose(result['LB'], expected_output["LB"])
    assert np.isclose(result['UB'], expected_output["UB"])

def test_rarity_MPN():
    MPN = 74.22876396932874
    positive = np.array([10, 5, 1])
    tubes = np.array([10, 10, 10])
    amount = np.array([0.1, 0.01, 0.001])

    expected_output = 0.7705316827363445

    result = rarity_MPN(MPN, positive, tubes, amount)

    assert isinstance(result, float)
    assert np.isclose(result, expected_output)

def test_isMissing():
    x1 = [1, 2, 3, float('nan'), 5]
    x2 = [1, 2, 3, float('inf'), 5]
    x3 = np.array([1, 2, 3, np.inf, 5])

    expected_output1 = True
    expected_output2 = True
    expected_output3 = True

    result1 = isMissing(x1)
    result2 = isMissing(x2)
    result3 = isMissing(x3)

    assert isinstance(result1, bool)
    assert isinstance(result2, bool)
    assert isinstance(result3, bool)
    assert result1 == expected_output1
    assert result2 == expected_output2
    assert result3 == expected_output3

def test_checkInputs_mpn():
    positive = [1, 3, 5, 0, 2]
    tubes = [5, 5, 5, 5, 5]
    amount = [1, 0.1, 0.01, 0.001, 0.0001]
    conf_level = 0.95

    # No exception should be raised for valid inputs
    checkInputs_mpn(positive, tubes, amount, conf_level)

    # Test cases for raising exceptions
    with pytest.raises(ValueError, match="positive, tubes, and amount must be the same length"):
        checkInputs_mpn(positive + [1], tubes, amount, conf_level)
    with pytest.raises(ValueError, match="'positive', 'tubes', and 'amount' cannot have missing values"):
        checkInputs_mpn([1, 3, 5, 0, float('nan')], tubes, amount, conf_level)
    with pytest.raises(TypeError, match="must be real number, not str"):
        checkInputs_mpn(positive, tubes, [1, 0.1, 0.01, 0.001, "abc"], conf_level)
    with pytest.raises(ValueError, match="'tubes' must contain positive whole numbers"):
        checkInputs_mpn(positive, [-1, 3, 5, 0, 2], amount, conf_level)
    with pytest.raises(ValueError, match="'positive' must contain non-negative whole numbers"):
        checkInputs_mpn([-1, 3, 5, 0, 2], tubes, amount, conf_level)
    with pytest.raises(ValueError, match="'amount' must contain positive values"):
        checkInputs_mpn(positive, tubes, [1, 0.1, -0.01, 0.001, 0.0001], conf_level)
    with pytest.raises(ValueError, match="'amount' must be in descending order"):
        checkInputs_mpn(positive, tubes, [1, 0.1, 0.01, 0.1, 0.0001], conf_level)
    with pytest.raises(ValueError, match="'conf_level' must be a float"):
        checkInputs_mpn(positive, tubes, amount, "0.95")
    with pytest.raises(ValueError, match="'conf_level' must be between 0 & 1"):
        checkInputs_mpn(positive, tubes, amount, -0.5)
        checkInputs_mpn(positive, tubes, amount, 1.5)
    with pytest.raises(ValueError, match="more positive tubes than possible"):
        checkInputs_mpn([1, 6, 5, 0, 2], tubes, amount, conf_level)


def test_mpn():
    positive = [10, 5, 1]
    tubes = [10, 10, 10]
    amount = [0.1, 0.01, 0.001]
    conf_level = 0.95
    CI_method = "Jarvis"

    expected_output = {'MPN': 74.22876396932874, 'MPN_adj': 67.92847966627644, 'variance': 902.4173735493894, 'var_log': 0.1637808007299598, 'conf_level': 0.95, 'CI_method': 'Jarvis', 'LB': 33.58089913036114, 'UB': 164.07867398144504, 'RI': 0.7705316827363445}


    result = mpn(positive, tubes, amount, conf_level, CI_method)

    assert isinstance(result, dict)
    assert len(result) == len(expected_output)
    assert np.isclose(result["MPN"], expected_output["MPN"])
    assert np.isclose(result["MPN_adj"], expected_output["MPN_adj"])
    assert np.isclose(result["variance"], expected_output["variance"])
    assert np.isclose(result["var_log"], expected_output["var_log"])
    assert np.isclose(result["conf_level"], expected_output["conf_level"])
    assert result["CI_method"] == expected_output["CI_method"]
    assert np.isclose(result["LB"], expected_output["LB"])
    assert np.isclose(result["UB"], expected_output["UB"])
    assert np.isclose(result["RI"], expected_output["RI"])

    CI_method = "LR"
    expected_output =  {'MPN': 74.22876396932874, 'MPN_adj': 67.92847966627644, 'variance': None, 'var_log': None, 'conf_level': 0.95, 'CI_method': 'LR', 'LB': 32.43393676359959, 'UB': 152.3661892483937, 'RI': 0.7705316827363445}

    result = mpn(positive, tubes, amount, conf_level, CI_method)

    assert isinstance(result, dict)
    assert len(result) == len(expected_output)
    assert np.isclose(result["MPN"], expected_output["MPN"])
    assert np.isclose(result["MPN_adj"], expected_output["MPN_adj"])
    assert result["variance"] is None
    assert result["var_log"] is None
    assert np.isclose(result["conf_level"], expected_output["conf_level"])
    assert result["CI_method"] == expected_output["CI_method"]
    assert np.isclose(result["LB"], expected_output["LB"])
    assert np.isclose(result["UB"], expected_output["UB"])
    assert np.isclose(result["RI"], expected_output["RI"])