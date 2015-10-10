#! /usr/bin/env python

"""
CLI program for estimating the posterior probability of divergence-time
scenarios from posterior samples of trees.
"""

import sys
import os
import re
import io
import collections
import gzip
import tempfile
import argparse
import logging
try:
    import queue
except ImportError:
    import Queue as queue
import multiprocessing

logging.basicConfig(level=logging.INFO)
_LOG = logging.getLogger(os.path.basename(__file__))
_LOCK = multiprocessing.Lock()

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1',
    'description': __doc__,
    'copyright': 'Copyright (C) 2015 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}


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


def is_gzipped(file_path):
    """
    Returns True if argument is a path to a gzipped file; False otherwise.

    Returns False if path does not exist:
    >>> is_gzipped('')
    False

    Returns False if regular file is empty:
    >>> fd, temp_path = tempfile.mkstemp()
    >>> os.close(fd)
    >>> is_gzipped(temp_path)
    False

    Returns False for regular file with content:
    >>> f = io.open(temp_path, mode = 'w', encoding='utf-8')
    >>> f.write('testing...')
    10
    >>> f.close()
    >>> is_gzipped(temp_path)
    False
    >>> os.remove(temp_path)

    Returns False if gzipped file is empty:
    >>> fd, temp_path = tempfile.mkstemp()
    >>> f = gzip.open(temp_path, mode = "wb", compresslevel = 9)
    >>> f.write(bytes("", 'UTF-8'))
    0
    >>> f.close()
    >>> is_gzipped(temp_path)
    False

    Returns True if file has gzipped content:
    >>> f = gzip.open(temp_path, mode = "wb", compresslevel = 9)
    >>> f.write(bytes("testing...", 'UTF-8'))
    10
    >>> f.close()
    >>> is_gzipped(temp_path)
    True
    >>> os.remove(temp_path)
    """

    try:
        fs = gzip.open(file_path)
        d = next(fs)
        fs.close()
    except:
        return False
    return True


class ReadFile(object):
    """
    Obtain a text stream in read mode from a regular or gzipped file.

    Behaves like ``open`` for regular files:
    >>> test_path = os.path.join(os.path.dirname(__file__), os.path.pardir,
    ...         'test-data', 'config.yml')
    >>> if os.path.exists(test_path):
    ...     with ReadFile(test_path) as f:
    ...         l = f.next().strip()
    ...     l == "---"
    ... else:
    ...     True
    ...
    ...
    True
    
    Behaves like ``open`` for gzipped files:
    >>> test_path = os.path.join(os.path.dirname(__file__), os.path.pardir,
    ...         'test-data', 'trees', 'crocs-1.trees.gz')
    >>> if os.path.exists(test_path):
    ...     with ReadFile(test_path) as f:
    ...         l = f.next().strip()
    ...     l == "#NEXUS"
    ... else:
    ...     True
    ...
    ...
    True
    """

    open_files = set()

    def __init__(self, path):
        self.path = path
        self.gzipped = is_gzipped(self.path)
        self.encoding = 'utf-8'
        if self.gzipped:
            self.file_stream = io.TextIOWrapper(
                    buffer = gzip.GzipFile(filename = self.path,
                            mode = 'rb'),
                    encoding = self.encoding)
        else:
            self.file_stream = io.open(self.path, mode = 'r',
                    encoding = self.encoding)
        self.__class__.open_files.add(self.path)

    def close(self):
        self.file_stream.close()
        self.__class__.open_files.remove(self.name)

    def __enter__(self):
        return self.file_stream

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        return self.file_stream.__iter__()

    def next(self):
        return next(self.file_stream)

    def read(self):
        return self.file_stream.read()

    def readline(self):
        return self.file_stream.readline()

    def readlines(self):
        return self.file_stream.readlines()

    def _get_closed(self):
        return self.file_stream.closed

    closed = property(_get_closed)

    def seek(self, i):
        self.file_stream.seek(i)

    def flush(self):
        self.file_stream.flush()

    def fileno(self):
        return self.file_stream.fileno()


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


class SumDivTimesError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class TipSubsetDataError(SumDivTimesError):
    def __init__(self, *args, **kwargs):
        SumDivTimesError.__init__(self, *args, **kwargs)


class PosteriorSampleDataError(SumDivTimesError):
    def __init__(self, *args, **kwargs):
        SumDivTimesError.__init__(self, *args, **kwargs)


class PosteriorManager(object):
    """
    Contains all PosteriorSamples, runs their workers and gets results
    back where they need to go
    """
    pass


class PosteriorSample(object):
    """
    This is the main container for the samples of divergence times from a
    posterior of phylogenies.

    An instance of this class should be initiated with parameters parsed from a
    ``posterior`` within the YAML config file. The parameters must include the
    following key-value pairs:

    -   paths : A list of strings that specify paths to tree files.
    -   schema : A string specifying the format of the tree files (default:
        'nexus').
    -   tip_subsets : A list of sets of key-value pairs that will each initiate
        a ``TipSubset`` instance.

    and optionally:

    -   name : A string.
    -   burnin : An integer.

    An example of initiating an instance:
    >>> d = {'paths': ['test-data/trees/test.trees.gz'],
    ...      'tip_subsets': [{
    ...             'name': 'bufo',
    ...             'tips': ['Bnebulifer', 'Bamericanus']}]}
    >>> ps = PosteriorSample(**d)
    >>> len(ps.paths) == 1
    True
    >>> ps.paths == (d['paths'][0],)
    True
    >>> len(ps.tip_subsets) == 1
    True
    >>> ps.tip_subsets[0].name == 'bufo'
    True
    """

    count = 0
    def __init__(self, *args, **kwargs):
        self.__class__.count += 1
        self.name = kwargs.pop('name',
                self.__class__.__name__ + '-' + str(self.count))
        paths = kwargs.pop('paths', None)
        if not paths:
            raise PosteriorSampleDataError("A posterior sample must contain a "
                    "list of 'paths'")
        self.paths = tuple(p for p in paths)
        self.schema = kwargs.pop('schema', None)
        if not schema:
            raise PosteriorSampleDataError("A posterior sample must specify "
                    "the 'schema' of the tree files")
        tip_subsets = kwargs.pop('tip_subsets', None)
        if not tip_subsets:
            raise PosteriorSampleDataError("A posterior sample must contain a "
                    "list of 'tip_subsets'")
        tip_subset_list = [TipSubset(**d) for d in tip_subsets]
        self.tip_subset_map = collections.OrderedDict(zip(
                [id(ts) for ts in tip_subset_list],
                tip_subset_list))
        self.burnin = int(kwargs.pop('burnin', 0))
        if len(kwargs) > 0:
            _LOG.warning("Unexpected attributes in posterior sample {0!r}: "
                    "{1}".format(self.name, ", ".join(kwargs.keys())))

    def _get_tip_subsets(self):
        return self.tip_subset_map.values()

    tip_subsets = property(_get_tip_subsets)


class PosteriorWorker(object):
    """
    These instances are assigned a single posterior tree file from which to
    extract divergence times
    """

    count = 0
    def __init__(self, path, schema, burnin, tip_subsets, label = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.path = path
        self.schema = schema
        self.burnin = burnin
        self.tip_subsets = tip_subsets
        self.label = label
        self.finished = False
        self.div_times = dict(zip(
                [id(ts) for ts in self.tip_subsets],
                [[] for i in range(len(self.tip_subsets))]))
        self.error = None
        self.trace_back = None

    def extract_div_times(self):
        try:
            return self._extract_div_times()
        except Exception as e:
            self.error = e
            f = StringIO()
            traceback.print_exc(file=f)
            self.trace_back = f.getvalue()


    def _extract_div_times(self):
        with ReadFile(self.path) as tree_stream:
            kwargs = {}
            if self.schema.lower() == 'nexus':
                kwargs['preserve_underscores'] = True
            tree_iter = dendropy.Tree.yield_from_files(
                    files = [tree_stream],
                    schema = self.schema,
                    tree_offset = self.burnin,
                    preserve_underscores = True)
            for tree_idx, tree in enumerate(tree_iter):
                if not tree.is_rooted:
                    raise Exception('Tree {0} in {1!r} is not rooted'.format(
                            i + 1,
                            self.path))
                for tip_subset in self.tip_subsets:
                    mrca_node = tree.mrca(taxon_labels = tip_subset.tips)
                    target_node = mrca_node
                    if target_node.parent_node and tip_subset.stem_based:
                        target_node = mrca_node.parent_node
                    self.div_times[id(tip_subset)].append((
                            id(target_node),
                            target_node.age))

class JobManager(multiprocessing.Process):
    count = 0
    def __init__(self,
            work_queue = None,
            result_queue = None,
            get_timeout = 0.4,
            put_timeout = 0.2,
            log = None,
            lock = None):
        multiprocessing.Process.__init__(self)
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        if not work_queue:
            work_queue = multiprocessing.Queue()
        if not result_queue:
            result_queue = multiprocessing.Queue()
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.get_timeout = get_timeout
        self.put_timeout = put_timeout
        if not log:
            log = _LOG
        self.log = log
        if not lock:
            lock = _LOCK
        self.lock = lock
        self.killed = False

    def compose_msg(self, msg):
        return '{0} ({1}): {2}'.format(self.name, self.pid, msg)

    def send_msg(self, msg, method_str='info'):
        self.lock.acquire()
        try:
            getattr(self.log, method_str)(self.compose_msg(msg))
        finally:
            self.lock.release()

    def send_debug(self, msg):
        self.send_msg(msg, method_str='debug')

    def send_info(self, msg):
        self.send_msg(msg, method_str='info')

    def send_warning(self, msg):
        self.send_msg(msg, method_str='warning')

    def send_error(self, msg):
        self.send_msg(msg, method_str='error')

    def _get_worker(self):
        worker = None
        try:
            self.send_debug('getting worker')
            worker = self.work_queue.get(block=True, timeout=self.get_timeout)
            self.send_debug('received worker {0}'.format(
                    getattr(worker, 'name', 'nameless')))
            # without blocking processes were stopping when the queue
            # was not empty, and without timeout, the processes would
            # hang waiting for jobs.
        except queue.Empty:
            time.sleep(0.2)
            if not self.work_queue.empty():
                self.send_warning('raised queue.Empty, but queue is '
                        'not empty... trying again')
                return self._get_worker()
            else:
                self.send_info('work queue is empty')
        return worker

    def _put_worker(self, worker):
        try:
            self.send_debug('returning worker {0}'.format(
                    getattr(worker, 'name', 'nameless')))
            self.result_queue.put(worker, block=True, timeout=self.put_timeout)
            self.send_debug('worker {0} returned'.format(
                    getattr(worker, 'name', 'nameless')))
        except queue.Full, e:
            time.sleep(0.2)
            if not self.result_queue.full():
                self.send_warning('raised queue.Full, but queue is '
                        'not full... trying again')
                self._put_worker(worker)
            else:
                self.send_error('result queue is full... aborting')
            self.killed = True
            raise e

    def run(self):
        self.send_debug('starting run')
        while not self.killed:
            worker = self._get_worker()
            if worker is None:
                break
            self.send_info('starting worker {0}'.format(
                    getattr(worker, 'name', 'nameless')))
            worker.extract_div_times()
            self.send_info('worker {0} finished'.format(
                    getattr(worker, 'name', 'nameless')))
            self._put_worker(worker)
        if self.killed:
            self.send_error('job manager was killed!')
        self.send_debug('end run')


    @classmethod
    def run_workers(cls,
            workers,
            num_processors,
            get_timeout = 0.4,
            put_timeout = 0.2,
            queue_max = 500):
        work_queue = multiprocessing.Queue()
        result_queue = multiprocessing.Queue()
        finished = []
        for w_list in list_splitter(workers, queue_max, by_size = True):
            assert work_queue.empty()
            assert result_queue.empty()
            for w in w_list:
                work_queue.put(w)
            managers = []
            for i in range(num_processors):
                m = cls(work_queue = work_queue,
                        result_queue = result_queue,
                        get_timeout = get_timeout,
                        put_timeout = put_timeout)
                managers.append(m)
            for i in range(len(managers)):
                managers[i].start()
            for i in range(len(w_list)):
                _LOG.debug('JobManager.run_workers: getting result...')
                w_list[i] = result_queue.get()
                _LOG.debug('JobManager.run_workers: got result {0}'.format(
                        getattr(w_list[i], 'name', 'nameless')))
            for i in range(len(managers)):
                managers[i].join()
            for w in w_list:
                if getattr(w, "error", None):
                    _LOG.error('Worker {0} returned with an error:\n{1}'.format(
                            getattr(w, 'name', 'nameless'),
                            w.trace_back))
                    raise w.error
            assert work_queue.empty()
            assert result_queue.empty()
            finished.extend(w_list)
        return finished


class TipSubset(object):
    """
    A subset of tips on trees from a posterior sample.

    An instance of this class should be initiated with parameters parsed from
    the items of a ``tip_subsets`` list within the YAML config file. The
    parameters must include the following key-value pairs:
    
    -   name : A string.
    -   tips : A list of strings that represent tip labels.

    and optionally:

    -   stem_based : A boolean (default ``False``).

    An example of initiating a TipSubset object:
    >>> d = {'name': 'bufo', 'tips': ['Bnebulifer', 'Bamericanus']}
    >>> ts = TipSubset(**d)
    >>> ts.name == d['name']
    True
    >>> ts.tips == tuple(t for t in d['tips'])
    True
    >>> ts.stem_based
    False

    'stem_based' keyword can also be passed:
    >>> d = {'name': 'bufo', 'tips': ['Bnebulifer', 'Bamericanus'], 'stem_based': True}
    >>> ts = TipSubset(**d)
    >>> ts.name == d['name']
    True
    >>> ts.tips == tuple(t for t in d['tips'])
    True
    >>> ts.stem_based
    True

    If 'name' or 'tips' is missing a TipSubsetDataError is raised:
    >>> d = {'tips': ['Bnebulifer', 'Bamericanus']}
    >>> ts = TipSubset(**d) #doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    TipSubsetDataError: ...

    >>> d = {'name': 'bufo'}
    >>> ts = TipSubset(**d) #doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    TipSubsetDataError: ...

    If any extra keywords are passed, a warning is given:
    >>> d = {'name': 'bufo', 'tips': ['Bnebulifer', 'Bamericanus'], 'extra': True}
    >>> ts = TipSubset(**d) #doctest: +ELLIPSIS
    """

    def __init__(self, *args, **kwargs):
        self.name = kwargs.pop('name', None)
        if not self.name:
            raise TipSubsetDataError("A tip subset must contain a 'name' "
                    "attribute")
        try:
            tips = list(kwargs.pop('tips', None))
        except:
            raise TipSubsetDataError("A tip subset must contain a list of "
                    "'tips'")
        self.tips = tuple(t for t in tips)
        self.stem_based = kwargs.pop('stem_based', False)
        if len(kwargs) > 0:
            _LOG.warning("Unexpected attributes in tip subset {0!r}: "
                    "{1}".format(self.name, ", ".join(kwargs.keys())))

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
