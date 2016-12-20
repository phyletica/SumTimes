#! /usr/bin/env python

"""
CLI program for estimating Bayes factors from SumTimes results.
"""

import sys
import os
import math
import time
import argparse
import unittest

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


def calculate_bayes_factors(
        sumtimes_posterior_results_path,
        sumtimes_prior_results_path):
    with open(sumtimes_posterior_results_path, 'r') as yaml_stream:
        posterior_results = yaml.load(yaml_stream)
        assert len(posterior_results) > 0
    with open(sumtimes_prior_results_path, 'r') as yaml_stream:
        prior_results = yaml.load(yaml_stream)
        assert len(prior_results) > 0

    posterior_expressions = {}
    prior_expressions = {}
    expected_keys = ['estimated_posterior_probability', 'expression', 'number_of_samples', 'number_of_samples_passing']
    for result_dict in posterior_results:
        assert sorted(result_dict.keys()) == expected_keys
        posterior_expressions[result_dict['expression'].strip()] = result_dict
    for result_dict in prior_results:
        assert sorted(result_dict.keys()) == expected_keys
        prior_expressions[result_dict['expression'].strip()] = result_dict

    posterior_expression_set = set(posterior_expressions.keys())
    prior_expression_set = set(prior_expressions.keys())
    if len(posterior_expression_set) < len(posterior_expressions):
        sys.stderr.write("ERROR: Found duplicated expressions in {0!r}".format(
                sumtimes_posterior_results_path))
        sys.exit(1)
    if len(prior_expression_set) < len(prior_expressions):
        sys.stderr.write("ERROR: Found duplicated expressions in {0!r}".format(
                sumtimes_prior_results_path))
        sys.exit(1)

    common_expression_set = posterior_expression_set.intersection(
            prior_expression_set)

    if len(common_expression_set) < 1:
        sys.stderr.write("ERROR: No common expressions in {0!r} and {1!r}".format(
                sumtimes_posterior_results_path,
                sumtimes_prior_results_path))
        sys.exit(1)

    bayes_factors = {}
    for expression in common_expression_set:
        posterior_k = posterior_expressions[expression]['number_of_samples_passing']
        posterior_n = posterior_expressions[expression]['number_of_samples']
        prior_k = prior_expressions[expression]['number_of_samples_passing']
        prior_n = prior_expressions[expression]['number_of_samples']

        posterior_odds_inf = False
        posterior_odds_zero = False
        if (posterior_k == 0):
            posterior_odds_zero = True
            posterior_odds = (posterior_k + 1) / float((posterior_n - posterior_k) - 1)
        elif (posterior_n == posterior_k):
            posterior_odds_inf = True
            posterior_odds = (posterior_k - 1) / float((posterior_n - posterior_k) + 1)
        else:
            posterior_odds = posterior_k / float(posterior_n - posterior_k)

        prior_odds_inf = False
        prior_odds_zero = False
        if (prior_k == 0):
            prior_odds_zero = True
            prior_odds = (prior_k + 1) / float((prior_n - prior_k) - 1)
        elif (prior_n == prior_k):
            prior_odds_inf = True
            prior_odds = (prior_k - 1) / float((prior_n - prior_k) + 1)
        else:
            prior_odds = prior_k / float(prior_n - prior_k)

        bf = posterior_odds / prior_odds
        lnbf = 2 * math.log(bf)

        bf_str = "undefined"
        lnbf_str = "undefined"
        if prior_odds_zero:
            if posterior_odds_zero:
                pass
            else:
                assert (not prior_odds_inf)
                bf_str = "> {0}".format(bf)
                lnbf_str = "> {0}".format(lnbf)
        else:
            if posterior_odds_zero:
                bf_str = "< {0}".format(bf)
                lnbf_str = "< {0}".format(lnbf)
            elif prior_odds_inf:
                if posterior_odds_inf:
                    pass
                else:
                    bf_str = "< {0}".format(bf)
                    lnbf_str = "< {0}".format(lnbf)
            else:
                if posterior_odds_inf:
                    bf_str = "> {0}".format(bf)
                    lnbf_str = "> {0}".format(lnbf)
                else:
                    bf_str = "{0}".format(bf)
                    lnbf_str = "{0}".format(lnbf)

        sys.stdout.write("- expression: >\n")
        sys.stdout.write("    {0}\n".format(expression))
        sys.stdout.write("  bf: {0}\n".format(bf_str))
        sys.stdout.write("  2lnbf: {0}\n".format(lnbf_str))
        bayes_factors[expression] = {
                'bf': bf,
                'lnbf': lnbf,
                'bf_str': bf_str,
                'lnbf_str': lnbf_str,
                }
    return bayes_factors


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
    bf = calculate_bayes_factors(
            args.sumtimes_posterior_results_path,
            args.sumtimes_prior_results_path)



class SumTimesBfCalculatorTestCase(unittest.TestCase):
    SCRIPT_PATH = os.path.abspath(__file__)
    BASE_DIR = os.path.abspath(os.path.dirname(SCRIPT_PATH))
    TEST_DATA_DIR = os.path.join(BASE_DIR, "test-data")

    def data_path(self, filename=""):
        return os.path.join(self.TEST_DATA_DIR, filename)
    
    def script_path(self):
        return self.SCRIPT_PATH

class IsFileTestCase(SumTimesBfCalculatorTestCase):
    def setUp(self):
        self.path = self.data_path('config.yml')
        self.bogus_path = self.data_path("bogusdatafilename")
    
    def test_is_file(self):
        self.assertFalse(is_file(None))
        self.assertFalse(is_file(self.bogus_path))
        self.assertTrue(is_file(self.path))

class ArgIsFileTestCase(SumTimesBfCalculatorTestCase):
    def test_simple(self):
        expected_path = os.path.abspath(__file__)
        self.assertTrue(expected_path == arg_is_file(__file__))

    def test_error(self):
        self.assertRaises(argparse.ArgumentTypeError,
                arg_is_file,
                "/this/path/probably/is/not/on/anyones/system")

class CalculateBayeFactorsTestCase(SumTimesBfCalculatorTestCase):
    def setUp(self):
        self.posterior_path = self.data_path('results.yml')
        self.prior_path = self.data_path('results-no-data.yml')

    def test_bf_calcs(self):
        bfs = calculate_bayes_factors(
                self.posterior_path,
                self.prior_path)
        self.assertEqual(len(bfs), 12)

        n1 = 20
        k1 = 20
        n2 = 10
        k2 = 0
        expected_bf = ((k1-1) / float(n1-k1+1)) / ((k2+1) / float(n2-k2-1))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 6.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 6.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 6.0)']['bf_str'],
                "> {0}".format(expected_bf),
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 6.0)']['lnbf_str'],
                "> {0}".format(expected_lnbf),
                )

        n1 = 10
        k1 = 10
        n2 = 20
        k2 = 20
        expected_bf = ((k1-1) / float(n1-k1+1)) / ((k2-1) / float(n2-k2+1))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 7.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 7.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 7.0)']['bf_str'],
                "undefined",
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 7.0)']['lnbf_str'],
                "undefined",
                )

        n1 = 20
        k1 = 20
        n2 = 10
        k2 = 2
        expected_bf = ((k1-1) / float(n1-k1+1)) / ((k2) / float(n2-k2))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 8.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 8.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 8.0)']['bf_str'],
                "> {0}".format(expected_bf),
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 8.0)']['lnbf_str'],
                "> {0}".format(expected_lnbf),
                )

        n1 = 20
        k1 = 0
        n2 = 10
        k2 = 10
        expected_bf = ((k1+1) / float(n1-k1-1)) / ((k2-1) / float(n2-k2+1))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 9.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 9.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 9.0)']['bf_str'],
                "< {0}".format(expected_bf),
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 9.0)']['lnbf_str'],
                "< {0}".format(expected_lnbf),
                )

        n1 = 20
        k1 = 10
        n2 = 10
        k2 = 10
        expected_bf = ((k1) / float(n1-k1)) / ((k2-1) / float(n2-k2+1))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10.0)']['bf_str'],
                "< {0}".format(expected_bf),
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10.0)']['lnbf_str'],
                "< {0}".format(expected_lnbf),
                )

        n1 = 10
        k1 = 7
        n2 = 20
        k2 = 0
        expected_bf = ((k1) / float(n1-k1)) / ((k2+1) / float(n2-k2-1))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 11.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 11.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 11.0)']['bf_str'],
                "> {0}".format(expected_bf),
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 11.0)']['lnbf_str'],
                "> {0}".format(expected_lnbf),
                )

        n1 = 10
        k1 = 0
        n2 = 10
        k2 = 3
        expected_bf = ((k1+1) / float(n1-k1-1)) / ((k2) / float(n2-k2))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 12.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 12.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 12.0)']['bf_str'],
                "< {0}".format(expected_bf),
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 12.0)']['lnbf_str'],
                "< {0}".format(expected_lnbf),
                )

        n1 = 20
        k1 = 0
        n2 = 10
        k2 = 0
        expected_bf = ((k1+1) / float(n1-k1-1)) / ((k2+1) / float(n2-k2-1))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 13.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 13.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 13.0)']['bf_str'],
                "undefined",
                )
        self.assertEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 13.0)']['lnbf_str'],
                "undefined",
                )

        n1 = 7501
        k1 = 4809
        n2 = 7500
        k2 = 3327
        expected_bf = (k1 / float(n1-k1)) / (k2 / float(n2-k2))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 5.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 5.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 5.0)']['bf_str']),
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 5.0)']['lnbf_str']),
                expected_lnbf,
                places = 10)

        n1 = 7501
        k1 = 107
        n2 = 7500
        k2 = 144
        expected_bf = (k1 / float(n1-k1)) / (k2 / float(n2-k2))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0)']['bf_str']),
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0)']['lnbf_str']),
                expected_lnbf,
                places = 10)

        n1 = 7501
        k1 = 45
        n2 = 7500
        k2 = 8
        expected_bf = (k1 / float(n1-k1)) / (k2 / float(n2-k2))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 3.0) & codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0) & ({tihamicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica})) & ({arabicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica}))']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 3.0) & codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0) & ({tihamicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica})) & ({arabicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica}))']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 3.0) & codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0) & ({tihamicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica})) & ({arabicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica}))']['bf_str']),
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 3.0) & codiverged(nodes = [{dhufarensis}, {ehrenbergii}, {felixarabica}], window = 3.0) & ({tihamicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica})) & ({arabicus} > max({dhufarensis}, {ehrenbergii}, {felixarabica}))']['lnbf_str']),
                expected_lnbf,
                places = 10)

        n1 = 7501
        k1 = 7055
        n2 = 7500
        k2 = 1032
        expected_bf = (k1 / float(n1-k1)) / (k2 / float(n2-k2))
        expected_lnbf = 2 * math.log(expected_bf)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10-30)']['bf'],
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10-30)']['lnbf'],
                expected_lnbf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10-30)']['bf_str']),
                expected_bf,
                places = 10)
        self.assertAlmostEqual(
                float(bfs['codiverged(nodes = [{tihamicus}, {arabicus}], window = 10-30)']['lnbf_str']),
                expected_lnbf,
                places = 10)


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
