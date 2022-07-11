"""
Unit and regression test for the enzymemap package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import enzymemap


def test_enzymemap_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "enzymemap" in sys.modules
