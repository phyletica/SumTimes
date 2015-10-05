#! /usr/bin/env python

"""
CLI program for estimating the posterior probability of divergence-time
scenarios from posterior samples of trees.
"""

import sys
import os
import re
import argparse

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1',
    'description': __doc__,
    'copyright': 'Copyright (C) 2015 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

class ConditionEvaluator(object):
    """
    >>> ce = ConditionEvaluator(
    ...        expression_str = "x < y < z",
    ...        variable_keys = ['x', 'y', 'z'])
    >>> str(ce)
    'x < y < z'

    >>> ce.evaluate({'x': 1.0, 'y': 2.0, 'z': 3.0})
    True

    Extra keys are okay:
    >>> ce.evaluate({'x': 1.0, 'y': 2.0, 'z': 3.0, 'xx': 4.0})
    True

    But, all variables in the expression must be keys:
    >>> ce.evaluate({'xx': 1.0, 'y': 2.0, 'z': 3.0})
    Traceback (most recent call last):
        ...
    KeyError: 'x'

    >>> ce = ConditionEvaluator("x < y = z", ['x', 'y', 'z'])
    Traceback (most recent call last):
        ...
    SyntaxError: Expression contains invalid syntax '='

    >>> ce = ConditionEvaluator("x < (y == z", ['x', 'y', 'z'])
    >>> ce.evaluate({'x': 1.0, 'y': 3.0, 'z': 3.0}) #doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    SyntaxError: ...

    >>> ce = ConditionEvaluator("abs(x - y) < 1.0", ['x', 'y'])
    >>> ce.evaluate({'x': 1.1, 'y': 1.3})
    True
    """

    INVALID_PATTERNS = {
            '=': re.compile(r'[^=]=[^=]'),
            }

    def __init__(self, expression_str, variable_keys):
        for k, p in self.INVALID_PATTERNS.items():
            if p.search(expression_str):
                raise SyntaxError('Expression contains invalid syntax '
                        '{0!r}'.format(k))
        self.variable_keys = list(variable_keys)
        self.expression_str = expression_str
        exp = self.expression_str
        for k in self.variable_keys:
            if exp.find(k) > -1:
                exp = exp.replace(k, "d['{0}']".format(k))
            else:
                raise KeyError('Expression does not contain key {0!r}'.format(v))
        self.expression = exp

    def __str__(self):
        return self.expression_str

    def evaluate(self, d):
        return eval(self.expression)


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
    Returns expanded path if its a file; returns argparse error otherwise

    >>> expected_path = os.path.abspath(__file__)
    >>> expected_path == arg_is_file(__file__)
    True

    >>> arg_is_file("/this/path/probably/is/not/on/anyones/system") #doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    argparse.ArgumentTypeError: ...
    """

    try:
        if not is_file(path):
            raise
    except:
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)


def main_cli(argv = sys.argv):
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument('config_path',
            metavar = 'YAML-CONFIG-PATH',
            type = arg_is_file,
            help = ('Path to the YAML-formatted config file.'))
    parser.add_argument('--run-tests',
            action = 'store_true',
            help = 'Ignore all other options and run test suite.')

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    if args.run_tests:
        import doctest
        doctest.testmod(verbose = True)

if __name__ == "__main__":
    main_cli()

