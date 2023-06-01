"""MicroMPN a Python command line program for automating microbiology most probable 
   number (MPN) estimates in laboratory microplates

    This module of MicroMPN is a transpiled version of a portion of the MPN R package:

    Martine Ferguson and John Ihrie. (2019). MPN: Most Probable Number 
    and Other Microbial Enumeration Techniques. R package version 0.3.0. 
    https://CRAN.R-project.org/package=MPN
"""

import math
from typing import Union, Dict, List

import numpy as np
from scipy.optimize import root_scalar
from scipy.stats import binom, chi2, norm


def ptEst_MPN(positive: np.array, tubes: np.array, amount: np.array) -> float:
    """Estimate MPN from observations of growth at different dilutions by maximizing the likelihood.

     Detailed FDA protocols are `here <https://www.fda.gov/food/laboratory-methods-food/bam-appendix-2-most-probable-number-serial-dilutions>`


    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :return: _description_
    :rtype: float
    """


    def scoreFnc(lambda_, positive, tubes, amount):
        # score function: Find root to maximize likelihood.
        LHS = np.sum((positive * amount) / (1 - np.exp(-lambda_ * amount)))
        RHS = np.sum(tubes * amount)
        return LHS - RHS

    all_negative = np.sum(positive) == 0
    all_positive = np.array_equal(positive, tubes)

    if all_negative:
        return 0
    elif all_positive:
        return float('inf')
    else:
        res = root_scalar(scoreFnc, bracket=[1e-04, 1e+14],
                          args=(positive, tubes, amount),
                          method='brentq', rtol=1e-12, maxiter=10000)
        return res.root


def ptEstAdj_MPN(MPN: float, positive: np.array, tubes: np.array, amount: np.array) -> float:
    """Adjust for bias in MPN estimate.

    Method from:

    Haas CN (1989). "Estimation of microbial densities from dilution count experiments" Applied and Environmental Microbiology 55(8), 1934-1942.
    Salama IA, Koch GG, Tolley DH. (1978) "On the estimation of the most probable number in a serial dilution technique." Communications in Statistics - Theory and Methods, 7(13), 1267-1281.

    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :return: A bias adjusted estimate of MPN
    :rtype: float
    """

    if MPN == 0:
        return 0
    elif np.isinf(MPN):
        return np.nan
    else:
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            # np warning context is used because bias correction method produces overflow values of 
            # cosh and sinh for large MPN that don't appreciably effect the fianl result but are to large to compute
            amount2 = amount ** 2
            amount3 = amount ** 3
            lambda_V = MPN * amount
            one_or_more = 1 - np.exp(-lambda_V)
            zi = tubes * one_or_more
            cosh_term = np.cosh(lambda_V) - 1
            D = np.sum((amount2 * zi) / (2 * cosh_term))
            wi1 = amount2 / (2 * (one_or_more ** 2) * (D ** 3))
            wi2 = (amount3 * zi * np.sinh(lambda_V)) / (cosh_term ** 2)
            wi2[np.isnan(wi2)] = 0  # Inf/Inf is NaN (denominator goes to infinity more quickly)
            wi2 = np.sum(wi2)
            wi3 = amount3 / (one_or_more * cosh_term * (D ** 2))
            wi = wi1 * wi2 - wi3
            return MPN - 0.5 * np.sum(wi * zi * np.exp(-lambda_V))


def logL_MPN(lambda_: float, positive: np.array, tubes: np.array, amount: np.array) -> float:
    """Estimate the Log-likelihood of the MPN estimate

    :param lambda_: MPN
    :type lambda_: float
    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :return: Log-likelihood
    :rtype: float
    """
    # log-likelihood
    lambda_amount = lambda_ * amount
    first_term = positive * np.log(1 - np.exp(-lambda_amount))
    second_term = (tubes - positive) * lambda_amount
    return np.sum(first_term - second_term)


def logLR_MPN(lambda_: float, lambda_hat: float, positive: np.array, tubes: np.array, amount: np.array) -> float:
    """Estimate the Log-likelihood ratio of the MPN estimate

    :param lambda_: MPN
    :type lambda_: float
    :param lambda_hat: MPN estimate
    :type lambda_hat: float
    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :return: Log-likelihood ratio
    :rtype: float
    """
    # log-likelihood ratio
    logL_alt = logL_MPN(lambda_hat, positive, tubes, amount)
    logL_null = logL_MPN(lambda_, positive, tubes, amount)
    return 2 * (logL_alt - logL_null)


def logLRroot_MPN(lambda_: float, lambda_hat: float, positive: np.array, tubes: np.array, amount: np.array, crit_val: float) -> float:
    """Find roots to get Log ratio confidence limits

    :param lambda_: MPN
    :type lambda_: float
    :param lambda_hat: MPN estimate
    :type lambda_hat: float
    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :param crit_val: The critical value for estimating the confidence limits
    :type crit_val: float
    :return: _description_
    :rtype: float
    """
    # Find roots to get LR confidence limits
    log_like_ratio = logLR_MPN(lambda_, lambda_hat, positive, tubes, amount)
    return log_like_ratio - crit_val



def jarvisCI_MPN(MPN: float, positive: np.array, tubes: np.array, amount: np.array, conf_level: float) -> Dict[str, float]:
    """Estimate the confidence interval of the MPN estimate using the method in Jarvis et al. (2010)

    Jarvis B, Wilrich C, Wilrich P-T (2010). "Reconsideration of the derivation of Most Probable Numbers, their standard deviations, confidence bounds and rarity values." Journal of Applied Microbiology, 109, 1660-1667. https://doi.org/10.1111/j.1365-2672.2010.04792.x


    :param MPN: MPN estimate
    :type MPN: float
    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :param conf_level: A scalar value between zero and one for the confidence level. Typically 0.95 (i.e., a 95 percent confidence interval).
    :type conf_level: float
    :return: A dictionary with the following values:
        variance: The estimated variance of the MPN estimate. If all tubes are positive or all negative, variance will be NA.
        var_logMPN: The estimated variance of the natural log of the MPN estimate using the Delta Method (see Jarvis et al.). If all tubes are positive or all negative, var_log will be NA. 
        LB: The lower bound of the confidence interval.
        UB: The upper bound of the confidence interval.
    :rtype: Dict[str, float]
    """
    all_negative = np.sum(positive) == 0
    all_positive = np.array_equal(positive, tubes)
    var_logMPN = np.nan
    variance = np.nan
    sig_level = 1 - conf_level

    if not all_negative and not all_positive:
        exp_term = np.exp(-amount * MPN)
        numer = positive * amount ** 2 * exp_term
        denom = (1 - exp_term) ** 2
        variance = 1 / np.sum(numer / denom)

    if all_negative:
        LB = 0
        UB = np.log(1 / sig_level) / np.sum(amount * tubes)
    elif all_positive:
        def fLB(LB):
            return np.log(1 / sig_level) + np.sum(tubes * np.log(1 - np.exp(-amount * LB)))

        result = root_scalar(fLB, bracket=[1e-02, 1e+03], method='brentq', maxiter=1e+04)
        LB = result.root
        UB = np.inf
    else:
        crit_val = norm.ppf(sig_level / 2, loc=0, scale=1)  # asym. normal
        var_logMPN = variance / (MPN ** 2)  # delta method
        SE_log = np.sqrt(var_logMPN)
        ME_log = crit_val * SE_log
        UB = MPN * np.exp(-ME_log)
        LB = MPN * np.exp(ME_log)

    return {'variance': variance, 'var_logMPN': var_logMPN, 'LB': LB, 'UB': UB}


def likeRatioCI_MPN(MPN: float, positive: np.array, tubes: np.array, amount: np.array, conf_level: float)-> Dict[str, float]:
    """Estimation of the confidence interval by the likelihood ratio method.

    For details see:

    Ridout MS (1994). "A comparison of confidence interval methods for dilution series experiments." Biometrics, 50(1), 289-296.

    :param MPN: MPN estimate
    :type MPN: float
    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :param conf_level: A scalar value between zero and one for the confidence level. Typically 0.95 (i.e., a 95 percent confidence interval).
    :type conf_level: float
    :return: A dictionary containing:
        LB: The lower bound of the confidence interval.
        UB: The upper bound of the confidence interval.
    :rtype: Dict[str, float]
    """
    # Likelihood ratio confidence limits
    # See Ridout (1994) - "A Comparison of CI Methods..."
    all_negative = np.sum(positive) == 0
    all_positive = np.array_equal(positive, tubes)
    crit_val = chi2.ppf(conf_level, df=1, loc=0, scale=1)
    if all_negative or all_positive:
        jarvis_bounds = jarvisCI_MPN(MPN, positive, tubes, amount, conf_level)
        LB = jarvis_bounds['LB']
        UB = jarvis_bounds['UB']
    else:
        result1 = root_scalar(logLRroot_MPN, bracket=[1e-04, MPN], method='brentq', maxiter=10000,
                              args=(MPN, positive, tubes, amount, crit_val))
        result2 = root_scalar(logLRroot_MPN, bracket=[MPN, 5 * MPN], method='brentq', maxiter=10000,
                              args=(MPN, positive, tubes, amount, crit_val))
        LB = result1.root
        UB = result2.root

    return {'LB': LB, 'UB': UB}



def rarity_MPN(MPN: float, positive: np.array, tubes: np.array, amount: np.array) -> float:
    """estimate rarity index from MPN

    If the Rarity Index is less than 1e-04, the experimental results are highly improbable. The researcher may consider running the experiment again and/or changing the dilution levels.

    Jarvis B, Wilrich C, Wilrich P-T (2010). "Reconsideration of the derivation of Most Probable Numbers, their standard deviations, 
    confidence bounds and rarity values." Journal of Applied Microbiology, 109, 1660-1667. https://doi.org/10.1111/j.1365-2672.2010.04792.x


    :param MPN: MPN estimate
    :type MPN: float
    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: np.array
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: np.array
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: np.array
    :return: The rarity index
    :rtype: float
    """

    all_negative = sum(positive) == 0
    all_positive = (positive == tubes).all()

    if all_negative or all_positive:
        return 1
    else:
        probs = 1 - np.exp(-MPN * amount)
        positive_max = np.minimum(np.floor(probs * (tubes + 1)), tubes)
        most_prob_L = binom.pmf(positive_max, tubes, probs)
        actual_L = binom.pmf(positive, tubes, probs)
        rarity = np.prod(actual_L) / np.prod(most_prob_L)
        return rarity

def isMissing(x: Union[list, np.array]) -> bool:
    """Check if any values in a vector are missing or infinite

    :param x: a Vector of values
    :type x: Union[List, np.array]
    :return: true is any value is missing of infinite
    :rtype: bool
    """

    return any([math.isnan(i) or math.isinf(i) for i in x])

def checkInputs_mpn(positive: Union[list, np.array], tubes: Union[list,np.array], amount: Union[list, np.array], conf_level: float) -> None:
    """Validations of input data to the MPN function

    :param positive:  a vector of the count of positive wells at a specific dilution
    :type positive: Union[list, np.array]
    :param tubes:  A vector of the total number of wells at a dilution
    :type tubes: Union[list,np.array]
    :param amount:  A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: Union[list, np.array]
    :param conf_level: The confidence level for estimating upper and lower bounds of MPN
    :type conf_level: float
    :raises ValueError: Any validation failure raises a value error

    """

    l_positive = len(positive)
    l_tubes = len(tubes)
    l_amount = len(amount)
    if l_positive != l_tubes or l_tubes != l_amount:
        raise ValueError("'positive', 'tubes', & 'amount' must be the same length")
    for datavect in [positive, tubes, amount]:
        if isMissing(datavect):
            raise ValueError("'positive', 'tubes', and  'amount' cannot have missing values")
    for datavect in [positive, tubes, amount]:
        if not all(isinstance(i, (int, float)) for i in datavect):
            raise TypeError("'positive', 'tubes', 'amount', must be numeric")
    if any([t < 1 for t in tubes]) or any([t != round(t) for t in tubes]):
        raise ValueError("'tubes' must contain positive whole numbers")
    if any([p < 0 for p in positive]) or any([p != round(p) for p in positive]):
        raise ValueError("'positive' must contain non-negative whole numbers")
    if any([a <= 0 for a in amount]):
        raise ValueError("'amount' must contain positive values")
    if len(amount) > 1:
        if max([amount[i] - amount[i+1] for i in range(len(amount)-1)]) <= 0:
            raise ValueError("'amount' must be in descending order")
    if not isinstance(conf_level, float):
        raise ValueError("'conf_level' must be a float")
    if conf_level <= 0 or conf_level >= 1:
        raise ValueError("'conf_level' must be between 0 & 1")
    if any([p > t for p,t in zip(positive, tubes)]):
        raise ValueError("more positive tubes than possible")

def mpn(positive: Union[List[int], np.array], tubes: Union[List[int], np.array], amount: Union[List[float], np.array], conf_level: float = 0.95, CI_method: str = "Jarvis") -> dict:
    """mpn calculates the Most Probable Number (MPN) point estimate and confidence interval for microbial concentrations. Also calculates Blodgett's (2002, 2005, 2010) Rarity Index (RI).

    This python code was transpiled from:

    Martine Ferguson and John Ihrie. (2019). MPN: Most Probable Number and Other Microbial Enumeration Techniques. R package version 0.3.0. https://CRAN.R-project.org/package=MPN .
    
    As an example, assume we start with 3g of undiluted inoculum per tube, then use a 10-fold dilution for 2 dilutions. We now have amount = 3 * c(1, .1, .01).
    When all tubes are negative, the point estimate of MPN is zero (same approach as Jarvis et al.). The point estimate for the BAM tables "is listed as less than the lowest MPN for an outcome with at least one positive tube" (App.2).
    When all tubes are positive, the point estimate for MPN is Inf (same approach as Jarvis et al.) since no finite maximum likelihood estimate (MLE) exists. The BAM tables "list the MPN for this outcome as greater than the highest MPN for an outcome with at least one negative tube" (App.2). Here, the difference is probably trivial since the sample should be further diluted if all tubes test positive.
    The bias adjustment for the point estimate uses the method of Salama et al. (1978). Also see Haas (1989).
    Currently, confidence intervals can only be calculated using the Jarvis (2010) or likelihood ratio (LR) approach (Ridout, 1994). The BAM tables use an alternate approach.
    We slightly modified Jarvis' approach when all tubes are positive or all are negative; we use \alpha instead of \alpha / 2α/2 since these are one-sided intervals.
      The Ridout (1994) LR approach uses the same technique (with \alpha) for these two extreme cases.
    If the Rarity Index is less than 1e-04, the experimental results are highly improbable. The researcher may consider running the experiment again and/or changing the dilution levels.

    References:

    Bacteriological Analytical Manual, 8th Edition, Appendix 2, https://www.fda.gov/Food/FoodScienceResearch/LaboratoryMethods/ucm109656.htm

    Blodgett RJ (2002). "Measuring improbability of outcomes from a serial dilution test." Communications in Statistics: Theory and Methods, 31(12), 2209-2223. https://www.tandfonline.com/doi/abs/10.1081/STA-120017222

    Blodgett RJ (2005). "Serial dilution with a confirmation step." Food Microbiology, 22(6), 547-552. https://doi.org/10.1016/j.fm.2004.11.017

    Blodgett RJ (2010). "Does a serial dilution experiment's model agree with its outcome?" Model Assisted Statistics and Applications, 5(3), 209-215. https://doi.org/10.3233/MAS-2010-0157

    Haas CN (1989). "Estimation of microbial densities from dilution count experiments" Applied and Environmental Microbiology 55(8), 1934-1942.

    Haas CN, Rose JB, Gerba CP (2014). "Quantitative microbial risk assessment, Second Ed." John Wiley & Sons, Inc., ISBN 978-1-118-14529-6.

    Jarvis B, Wilrich C, Wilrich P-T (2010). "Reconsideration of the derivation of Most Probable Numbers, their standard deviations, confidence bounds and rarity values." Journal of Applied Microbiology, 109, 1660-1667. https://doi.org/10.1111/j.1365-2672.2010.04792.x

    Ridout MS (1994). "A comparison of confidence interval methods for dilution series experiments." Biometrics, 50(1), 289-296.

    Salama IA, Koch GG, Tolley DH. (1978) "On the estimation of the most probable number in a serial dilution technique." Communications in Statistics - Theory and Methods, 7(13), 1267-1281.



    :param positive: A vector of number of positive tubes at each dilution level.
    :type positive: Union[List, np.array]
    :param tubes: A vector of total number of tubes at each dilution level.
    :type tubes: Union[List,np.array]
    :param amount: A vector of the dilution for each group of wells (smaller is more dilute) 
    :type amount: Union[List, np.array]
    :param conf_level: A scalar value between zero and one for the confidence level. Typically 0.95 (i.e., a 95 percent confidence interval). Defaults to 0.95
    :type conf_level: float, optional
    :param CI_method: The method used for calculating the confidence interval. Choices are "Jarvis" or "LR" (likelihood ratio). Defaults to "Jarvis"
    :type CI_method: str, optional
    :return: A dictionary with these values:
             MPN: The most probable number point estimate for microbial density (concentration).
             MPN_adj: The bias-adjusted point estimate for MPN.
             variance: The estimated variance (see Jarvis et al.) of the MPN estimate if CI_method = "Jarvis". If all tubes are positive or all negative, variance will be NA. If CI_method is not "Jarvis", variance will be NA.
             var_log: The estimated variance of the natural log of the MPN estimate (see Jarvis et al.) using the Delta Method. If all tubes are positive or all negative, var_log will be NA. If CI_method is not "Jarvis", var_log will be NA.
             conf_level: The confidence level used.
             CI_method: The confidence interval method used.
             LB: The lower bound of the confidence interval.
             UB: The upper bound of the confidence interval.
             RI: The rarity index.
    :rtype: Dict[str, float]
    """
    checkInputs_mpn(positive=positive, tubes=tubes, amount=amount, conf_level=conf_level)
    positive_np = np.array(positive)
    tubes_np = np.array(tubes)
    amount_np = np.array(amount)
    MPN = ptEst_MPN(positive=positive_np, tubes=tubes_np, amount=amount_np)
    MPN_adj = ptEstAdj_MPN(MPN=MPN, positive=positive_np, tubes=tubes_np, amount=amount_np)

    variance = None
    var_logMPN = None
    if CI_method == "Jarvis":
        jarvis = jarvisCI_MPN(MPN=MPN, positive=positive_np, tubes=tubes_np, amount=amount_np, conf_level=conf_level)
        variance = jarvis["variance"]
        var_logMPN = jarvis["var_logMPN"]
        LB = jarvis["LB"]
        UB = jarvis["UB"]
    elif CI_method == "LR":
        like_ratio = likeRatioCI_MPN(MPN=MPN, positive=positive_np, tubes=tubes_np, amount=amount_np, conf_level=conf_level)
        LB = like_ratio["LB"]
        UB = like_ratio["UB"]

    rarity = rarity_MPN(MPN, positive_np, tubes_np, amount_np)

    return {"MPN": MPN, "MPN_adj": MPN_adj, "variance": variance, "var_log": var_logMPN,
            "conf_level": conf_level, "CI_method": CI_method, "LB": LB, "UB": UB, "RI": rarity}

