import pytest
import numpy as np
from src.analysis.residue_analysis import bfactor_regression
from src.analysis.validation import validate_bfactor_correlation

def test_bfactor():
    disp = np.random.randn(10, 5)
    bf = bfactor_regression(disp)
    assert len(bf) == 5
    assert np.all(bf > 0)

def test_validation():
    comp = np.arange(10)
    ref = comp * 0.8 + np.random.randn(10)*0.1
    assert validate_bfactor_correlation(comp, ref)
