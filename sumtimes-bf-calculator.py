#! /usr/bin/env python

"""
CLI program for estimating Bayes factors from SumTimes results.
"""

import sys
import os
import time
import argparse
import logging
import unittest

logging.basicConfig(level=logging.INFO)
_LOG = logging.getLogger(os.path.basename(__file__))
_LOCK = multiprocessing.Lock()
_RNG = random.Random()

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': '0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2015 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}


try:
    import yaml
except ImportError:
    raise ImportError("""No module named yaml.
The package 'pyyaml' is required for {0}.  Please try to install on the command
line using:

    sudo pip install pyyaml

If you do not have admin privileges you can install via:

    pip install --user pyyaml
""".format(_program_info['name']))


def expand_path(path):
    """
    Returns a full path

    >>> expected_path = os.path.expanduser('~/testing')
    >>> expected_path == expand_path('~/testing')
    True
    """

    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))


def is_file(path):
    """
    Returns boolean of whether or not argument is a path

    Returns False if path does not exist:
    >>> is_file("/this/path/probably/is/not/on/anyones/system")
    False

    Returns False if the path is a directory:
    >>> is_file(os.path.dirname(__file__))
    False
    
    Returns True if the path is a file:
    >>> is_file(__file__)
    True
    """

    if not path:
        return False
    if not os.path.isfile(path):
        return False
    return True

def arg_is_file(path):
    """
    Returns expanded path if argument is a file; returns argparse error
    otherwise.

    >>> expected_path = os.path.abspath(__file__)
    >>> expected_path == arg_is_file(__file__)
    True
    """

    try:
        if not is_file(path):
            raise
    except:
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)


def main_cli(argv = sys.argv):
    description = '{name} Version {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument('sumtimes_posterior_results_path',
            metavar = 'SUMTIMES-POSTERIOR-RESULTS-PATH',
            type = arg_is_file,
            help = ('Path to the SumTimes results path with estimated '
                    'posterior probabilities of divergence scenarios.'))
    parser.add_argument('sumtimes_prior_results_path',
            metavar = 'SUMTIMES-PRIOR-RESULTS-PATH',
            type = arg_is_file,
            help = ('Path to the SumTimes results path with estimated '
                    'prior probabilities of divergence scenarios.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)
    
    with open(args.sumtimes_posterior_results_path, 'r') as yaml_stream:
        posterior_results = yaml.load(yaml_stream)
        assert len(posterior_results) > 0
    with open(args.sumtimes_prior_results_path, 'r') as yaml_stream:
        prior_results = yaml.load(yaml_stream)
        assert len(prior_results) > 0

    posterior_expressions = []
    prior_expressions = []
    expected_keys = ['estimated_posterior_probability', 'expression', 'number_of_samples', 'number_of_samples_passing']
    for result_dict in posterior_results:
        assert sorted(result_dict.keys()) == expected_keys
        posterior_expressions.append(result_dict['expression'])
    for result_dict in prior_results:
        assert sorted(result_dict.keys()) == expected_keys
        prior_expressions.append(result_dict['expression'])


if __name__ == "__main__":
    if "--run-tests" in sys.argv:

        sys.stderr.write("""
*********************************************************************
Running test suite using the following Python executable and version:
{0}
{1}
*********************************************************************
\n""".format(sys.executable, sys.version))

        import doctest

        # doctest.testmod(verbose = True)
        suite = unittest.TestSuite()
        suite.addTest(doctest.DocTestSuite())

        tests = unittest.defaultTestLoader.loadTestsFromName(
               os.path.splitext(os.path.basename(__file__))[0])
        suite.addTests(tests)

        runner = unittest.TextTestRunner(verbosity = 2)
        runner.run(suite)

        sys.exit(0)

    main_cli()
