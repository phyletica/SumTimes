#! /usr/bin/env python

"""
CLI program for estimating the posterior probability of divergence-time
scenarios from posterior samples of trees.
"""

import sys
import os
import re
import io
import time
import collections
import gzip
import traceback
import tempfile
import argparse
import logging
import unittest
try:
    import queue
except ImportError:
    import Queue as queue
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import multiprocessing

logging.basicConfig(level=logging.INFO)
_LOG = logging.getLogger(os.path.basename(__file__))
_LOCK = multiprocessing.Lock()

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
    _LOG.error("""
Could not import yaml. The package 'pyyaml' is required for {0}.
Please try to install on the command line using:

    sudo pip install pyyaml

If you do not have admin privileges you can install via:

    pip install --user pyyaml
""".format(_program_info['name']))
    sys.exit(1)

try:
    import dendropy
except ImportError:
    _LOG.error("""
Could not import dendropy. The package 'dendropy' is required for {0}.
Please try to install on the command line using:

    sudo pip install dendropy

If you do not have admin privileges you can install via:

    pip install --user dendropy
""".format(_program_info['name']))
    sys.exit(1)


def list_splitter(l, n, by_size=False):
    """
    Returns generator that yields list `l` as `n` sublists, or as `n`-sized
    sublists if `by_size` is True.

    >>> import types
    >>> ret = list_splitter([1,2,3,4,5,6,7], 2)
    >>> isinstance(ret, types.GeneratorType)
    True
    >>> ret = list(ret)
    >>> len(ret) == 2
    True
    >>> isinstance(ret[0], list)
    True
    >>> isinstance(ret[1], list)
    True
    >>> ret[0] == [1,2,3] 
    True

    When ``by_size = False`` and length of input list is not multiple of ``n``,
    extra elements are in the last list.
    >>> ret[1] == [4,5,6,7]
    True

    >>> ret = list_splitter([1,2,3,4,5,6,7], 2, by_size = True)
    >>> isinstance(ret, types.GeneratorType)
    True
    >>> ret = list(ret)
    >>> len(ret) == 4
    True
    >>> isinstance(ret[0], list)
    True
    >>> ret[0] == [1,2] 
    True
    >>> ret[1] == [3,4] 
    True
    >>> ret[3] == [7] 
    True
    """
    if n < 1:
        raise StopIteration
    elif by_size:
        for i in range(0, len(l), n):
            yield l[i:i+n]
    else:
        if n > len(l):
            n = len(l)
        step_size = len(l)//int(n)
        if step_size < 1:
            step_size = 1
        i = -step_size
        for i in range(0, ((n-1)*step_size), step_size):
            yield l[i:i+step_size]
        yield l[i+step_size:]


def get_new_path(path, max_attempts = 1000):
    path = os.path.abspath(os.path.expandvars(os.path.expanduser(path)))
    if not os.path.exists(path):
        return path
    attempt = 0
    while True:
        p = '-'.join([path, str(attempt)])
        if not os.path.exists(p):
            return p
        if attempt >= max_attempts:
            raise Exception('failed to get unique path')
        attempt += 1


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
    >>> f.write(u'testing...') #doctest: +ELLIPSIS
    10...
    >>> f.close()
    >>> is_gzipped(temp_path)
    False
    >>> os.remove(temp_path)

    Returns False if gzipped file is empty:
    >>> fd, temp_path = tempfile.mkstemp()
    >>> f = gzip.open(temp_path, mode = "wt", compresslevel = 9)
    >>> f.write("")
    0
    >>> f.close()
    >>> is_gzipped(temp_path)
    False

    Returns True if file has gzipped content:
    >>> f = gzip.open(temp_path, mode = "wt", compresslevel = 9)
    >>> f.write("testing...")
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
            try:
                self.file_stream = gzip.open(filename = self.path, mode = 'rt',
                        encoding = self.encoding)
            except (TypeError, ValueError):
                self.file_stream = gzip.open(filename = self.path, mode = 'r')
        else:
            self.file_stream = io.open(self.path, mode = 'r',
                    encoding = self.encoding)
        self.__class__.open_files.add(self.path)

    def close(self):
        self.file_stream.close()
        self.__class__.open_files.remove(self.path)

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
        self.variable_keys = sorted(variable_keys, key = len)
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


class YamlConfigFormattingError(SumDivTimesError):
    def __init__(self, *args, **kwargs):
        SumDivTimesError.__init__(self, *args, **kwargs)


class AnalysisManager(object):
    """
    Container for all posterior samples found in the YAML config file.
    
    An instance is initiated with a list of ``posterior``s parsed from
    the YAML config file.
    """

    number_pattern_string = r'[\d\.Ee\-\+]+'
    codiverged_pattern_string = (
            '(?P<block>codiverged\s*\(\s*'
            'nodes\s*=\s*\[\s*(?P<nodes>.+)\s*\]'
            '\s*,\s*'
            'window\s*=\s*(?P<window>.+)\s*\))')
    codiverged_pattern = re.compile(codiverged_pattern_string)

    window_pattern_string = '(?P<min>{0})\s*-\s*(?P<max>{0})'.format(number_pattern_string)
    window_pattern = re.compile(window_pattern_string)

    def __init__(self, config_path, num_processors):
        posterior_samples = []
        self.config_path = config_path
        with ReadFile(self.config_path) as yaml_stream:
            config = yaml.load(yaml_stream)
        expression_strings = []
        for setting_dict in config:
            if len(setting_dict) > 1:
                raise YamlConfigFormattingError(
                        "Top level keys must be either 'posterior' or "
                        "'expression'; found: {0} ".format(", ".join(
                                setting_dict.keys())))
            key = list(setting_dict.keys())[0]
            if key.lower() == 'posterior':
                posterior_samples.append(PosteriorSample(
                        config_path = config_path,
                        **setting_dict[key]))
            elif key.lower() == 'expression':
                expression_strings.append(setting_dict[key])
            else:
                raise YamlConfigFormattingError(
                        "Top level keys must be either 'posterior' or "
                        "'expression'; found: {0} ".format(key))
        self.posterior_sample_map = collections.OrderedDict(zip(
                [ps.name for ps in posterior_samples],
                posterior_samples))
        self.num_processors = int(num_processors)
        self.shared_nodes = []
        tip_subset_names = []
        for posterior_sample in self.posterior_sample_map.values():
            for tip_subset in posterior_sample.tip_subset_map.values():
                tip_subset_names.append(tip_subset.name)
        self.tip_subset_names = set(tip_subset_names)
        for expression_str in expression_strings:
            self._parse_expression(expression_str)
        # self.shared_node_path = os.path.join(
        #         os.path.dirname(self.config_path),
        #         'shared-node-info.txt')


    # codiverged(nodes = [crocodylus, paleosuchus, mindorensis], window=8-12)
    # & (crocodylus > alligator)

    def _get_posterior_samples(self):
        return self.posterior_sample_map.values()

    posterior_samples = property(_get_posterior_samples)

    def _parse_expression(self, expression):
        for codiv_match in self.codiverged_pattern.finditer(expression):
            nodes_str = codiv_match.group('nodes')
            node_list = [n.strip() for n in nodes_str.split(',')]
            nodes = set(node_list)
            if len(nodes) != len(node_list):
                raise YamlConfigFormattingError(
                        "duplicated tip subset names in codiverged "
                        "expression:\n{0}".format(
                                codiv_match.group('block')))
            if not nodes.issubset(self.tip_subset_names):
                raise YamlConfigFormattingError(
                        "undefined tip subset names in codiverged "
                        "expression:\n{0}".format(
                                codiv_match.group('block')))
                
            window_str = codiv_match.group('window')
            window_width = None
            window_min = None
            window_max = None
            try:
                window_width = float(window_str)
            except:
                window_match = self.window_pattern.match(window_str)
                try:
                    window_min = float(window_match.group('min'))
                    window_max = float(window_match.group('max'))
                except:
                    raise YamlConfigFormattingError("could not parse codiverged"
                            " 'window' argument in:\n{0}".format(
                                    codiv_match.group('block')))
            replacement = StringIO()
            if window_width is not None:
                pass
            else:
                pass
            pass

    def _get_node_age_extractors(self):
        workers = []
        for ps in self.posterior_samples:
            workers.extend(ps.get_workers())
        return workers

    def _extract_node_ages(self):
        workers = self._get_node_age_extractors()
        np = min((self.num_processors, len(workers)))
        workers = JobManager.run_workers(workers, np)
        for worker in workers:
            for tip_subset_name, node_ages in worker.node_ages.items():
                self.posterior_sample_map[worker.label].tip_subset_map[
                        tip_subset_name].node_ages.extend(node_ages)
            self.shared_nodes.extend(worker.shared_nodes)

    def run_analysis(self):
        self._extract_node_ages()
        self._write_shared_node_age_warning()

    def _write_shared_node_age_warning(self):
        if self.shared_nodes:
            out_stream = StringIO()
            sys.stderr.write("""\
WARNING: Some of the tip subsets you defined map to the same node in some of
         the trees. All of the information about the shared nodes is listed
         below. For each case in which multiple tip subsets mapped to the same
         node, the tree file, index, and label is given, along with the tip
         subsets that shared a node.\n""")
            out_stream.write("tree_file\ttree_index\ttree_label\ttip_subsets\n")
            for shared_node in self.shared_nodes:
                out_stream.write("{0.tree_path}\t{0.tree_index}\t"
                        "{0.tree_label}\t{1}\n".format(
                                shared_node, 
                                ", ".join(shared_node.tip_subset_names)))
            sys.stderr.write(out_stream.getvalue())


class PosteriorSample(object):
    """
    This is the main container for the samples of node ages from a
    posterior of phylogenies.

    An instance of this class should be initiated with parameters parsed from a
    ``posterior`` within the YAML config file. One positional argument is also
    required:

    -   config_path : A string of the path to the YAML config file.
    
    The parameters must include the following key-value pairs:

    -   paths : A list of strings that specify paths to tree files (relative to
        the config path).
    -   tip_subsets : A list of sets of key-value pairs that will each initiate
        a ``TipSubset`` instance.

    and optionally:

    -   burnin : An integer.
    -   schema : A string specifying the format of the tree files (default:
        'nexus').

    An example of initiating an instance:
    >>> d = {'paths': ['trees/test.trees.gz'],
    ...      'tip_subsets': [{
    ...             'name': 'bufo',
    ...             'tips': ['Bnebulifer', 'Bamericanus']}]}
    >>> ps = PosteriorSample('test-data/config.yml', **d)
    >>> len(ps.paths) == 1
    True
    >>> ps.paths == (os.path.abspath('test-data/trees/test.trees.gz'),)
    True
    >>> len(ps.tip_subsets) == 1
    True
    >>> list(ps.tip_subsets)[0].name == 'bufo'
    True
    """

    count = 0

    def __init__(self, config_path, **kwargs):
        self.__class__.count += 1
        self.config_path = config_path
        self.name = self.__class__.__name__ + '-' + str(self.count)
        paths = kwargs.pop('paths', None)
        if not paths:
            raise PosteriorSampleDataError("A posterior sample must contain a "
                    "list of 'paths'")
        self.paths = tuple(expand_path(os.path.join(os.path.dirname(
                config_path), p)) for p in paths)
        tip_subsets = kwargs.pop('tip_subsets', None)
        if not tip_subsets:
            raise PosteriorSampleDataError("A posterior sample must contain a "
                    "list of 'tip_subsets'")
        tip_subset_list = [TipSubset(**d) for d in tip_subsets]
        self.tip_subset_map = collections.OrderedDict(zip(
                [ts.name for ts in tip_subset_list],
                tip_subset_list))
        self.schema = kwargs.pop('schema', 'nexus')
        self.burnin = int(kwargs.pop('burnin', 0))
        if len(kwargs) > 0:
            _LOG.warning("Unexpected attributes in posterior sample {0!r}: "
                    "{1}".format(self.name, ", ".join(kwargs.keys())))

    def _get_tip_subsets(self):
        return self.tip_subset_map.values()

    tip_subsets = property(_get_tip_subsets)

    def get_workers(self):
        workers = []
        for path in self.paths:
            w = PosteriorWorker(path = path,
                    schema = self.schema,
                    burnin = self.burnin,
                    tip_subsets = list(self.tip_subsets),
                    label = self.name)
            workers.append(w)
        return workers


class PosteriorWorker(object):
    """
    These instances are assigned a single posterior tree file from which to
    extract node ages.

    Instances should be created via the ``get_workers`` method of the
    ``PosteriorSample`` class

    >>> d = {'paths': ['trees/test.trees.gz'],
    ...      'tip_subsets': [{
    ...             'name': 'toads',
    ...             'tips': ['Bnebulifer', 'Bamericanus']}]}
    >>> ps = PosteriorSample('test-data/config.yml', **d)
    >>> workers = ps.get_workers()
    >>> len(workers) == 1
    True
    >>> w = workers[0]
    >>> w.path == os.path.abspath('test-data/trees/test.trees.gz')
    True
    >>> w.burnin == 0
    True
    >>> w.schema == 'nexus'
    True
    >>> w.label == ps.name
    True
    >>> w.finished == False
    True
    >>> w.node_ages == {'toads': []}
    True
    >>> w.tip_subsets == [('toads', ('Bnebulifer', 'Bamericanus'), False)]
    True
    """

    count = 0
    def __init__(self, path, schema, burnin, tip_subsets, label = None):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.path = path
        self.schema = schema
        self.burnin = burnin
        self.label = label
        self.finished = False
        self.node_ages = dict((ts.name, []) for ts in tip_subsets)
        self.tip_subsets = [(ts.name, ts.tips, ts.stem_based) for ts in tip_subsets]
        self.shared_nodes = []
        self.error = None
        self.trace_back = None

    def start(self):
        try:
            return self._extract_node_ages()
        except Exception as e:
            self.error = e
            f = StringIO()
            traceback.print_exc(file=f)
            self.trace_back = f.getvalue()


    def _extract_node_ages(self):
        with ReadFile(self.path) as tree_stream:
            tree_iter = dendropy.Tree.yield_from_files(
                    files = [tree_stream],
                    schema = self.schema,
                    preserve_underscores = True)
            for tree_idx, tree in enumerate(tree_iter):
                if tree_idx < self.burnin:
                    continue
                if not tree.is_rooted:
                    raise Exception('Tree {0} in {1!r} is not rooted'.format(
                            i + 1,
                            self.path))
                tree.calc_node_ages()
                nodes_visited = {}
                shared_nodes = set()
                for name, tips, stem_based in self.tip_subsets:
                    mrca_node = tree.mrca(taxon_labels = tips)
                    target_node = mrca_node
                    if target_node.parent_node and stem_based:
                        target_node = mrca_node.parent_node
                    self.node_ages[name].append(target_node.age)
                    if id(target_node) in nodes_visited:
                        nodes_visited[id(target_node)].append(name)
                        shared_nodes.add(id(target_node))
                    else:
                        nodes_visited[id(target_node)] = [name]
                for node_id in shared_nodes:
                    self.shared_nodes.append(SharedNode(
                            tree_path = self.path,
                            tree_index = tree_idx,
                            tree_label = tree.label,
                            node_id = node_id,
                            tip_subset_names = nodes_visited[node_id]))
        self.finished = True

class SharedNode(object):
    def __init__(self,
            tree_path,
            tree_index,
            tree_label,
            node_id,
            tip_subset_names):
        self.tree_path = tree_path
        self.tree_index = tree_index
        self.tree_label = tree_label
        self.node_id = node_id
        self.tip_subset_names = tuple(tip_subset_names)

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
        except queue.Full as e:
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
            worker.start()
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
    >>> d = {'name': 'buf', 'tips': ['Bnebulifer', 'Bamericanus']}
    >>> ts = TipSubset(**d)
    >>> ts.name == d['name']
    True
    >>> ts.tips == tuple(t for t in d['tips'])
    True
    >>> ts.stem_based
    False

    'stem_based' keyword can also be passed:
    >>> d = {'name': 'buf2', 'tips': ['Bnebulifer', 'Bamericanus'], 'stem_based': True}
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

    Every instance must have a unique name:
    >>> d = {'name': 'buf', 'tips': ['Bnebulifer', 'Bamericanus'], 'extra': True}
    >>> ts = TipSubset(**d) #doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    TipSubsetDataError: All tip subset names must be unique

    If any extra keywords are passed, a warning is given:
    >>> d = {'name': 'buf3', 'tips': ['Bnebulifer', 'Bamericanus'], 'extra': True}
    >>> ts = TipSubset(**d) #doctest: +ELLIPSIS
    """
    
    count = 0
    registered_names = set()

    def __init__(self, *args, **kwargs):
        self.__class__.count += 1
        name = kwargs.pop('name', None)
        if not name:
            raise TipSubsetDataError("A tip subset must contain a 'name' "
                    "attribute")
        if name in self.registered_names:
            raise TipSubsetDataError("All tip subset names must be unique")
        self.__class__.registered_names.add(name)
        self.name = name
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
        self.node_ages = []

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

def arg_is_positive_int(i):
    """
    Returns int if argument is a positive integer; returns argparse error
    otherwise.

    >>> arg_is_positive_int(1) == 1
    True
    """

    try:
        if int(i) < 1:
            raise
    except:
        msg = '{0!r} is not a positive integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)


def main_cli(argv = sys.argv):
    description = '{name} Version {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument('config_path',
            metavar = 'YAML-CONFIG-PATH',
            type = arg_is_file,
            help = ('Path to the YAML-formatted config file.'))
    parser.add_argument('--np',
            action = 'store',
            type = arg_is_positive_int,
            default = multiprocessing.cpu_count(),
            help = ('The maximum number of processes to run in parallel. The '
                    'default is the smaller of the number of tree files or the '
                    'number of CPUs available on the machine.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    analysis = AnalysisManager(config_path = args.config_path,
            num_processors = args.np)
    analysis.run_analysis()

    
class SumDivTimesTestCase(unittest.TestCase):
    SCRIPT_PATH = os.path.abspath(__file__)
    BASE_DIR = os.path.abspath(os.path.dirname(SCRIPT_PATH))
    TEST_DATA_DIR = os.path.join(BASE_DIR, "test-data")

    def data_path(self, filename=""):
        return os.path.join(self.TEST_DATA_DIR, filename)
    
    def script_path(self):
        return self.SCRIPT_PATH

class IsFileTestCase(SumDivTimesTestCase):
    def setUp(self):
        self.path = self.data_path('config.yml')
        self.bogus_path = self.data_path("bogusdatafilename")
    
    def test_is_file(self):
        self.assertFalse(is_file(None))
        self.assertFalse(is_file(self.bogus_path))
        self.assertTrue(is_file(self.path))

class ArgIsFileTestCase(SumDivTimesTestCase):
    def test_simple(self):
        expected_path = os.path.abspath(__file__)
        self.assertTrue(expected_path == arg_is_file(__file__))

    def test_error(self):
        self.assertRaises(argparse.ArgumentTypeError,
                arg_is_file,
                "/this/path/probably/is/not/on/anyones/system")

class ArgIsPositiveIntTestCase(SumDivTimesTestCase):
    def test_simple(self):
        self.assertEqual(arg_is_positive_int(1), 1)

    def test_error(self):
        self.assertRaises(argparse.ArgumentTypeError,
                arg_is_positive_int,
                -1)

class AnalysisManagerTestCase(SumDivTimesTestCase):
    def setUp(self):
        self.config_path = self.data_path('config.yml')
        TipSubset.registered_names = set()
        TipSubset.count = 0
        PosteriorSample.count = 0

    def test_init_simple(self):
        analysis = AnalysisManager(
                config_path = self.config_path,
                num_processors = 8)
        self.assertEqual(analysis.num_processors, 8)
        self.assertEqual(len(analysis.posterior_samples), 2)
        self.assertEqual(analysis.shared_nodes, [])

        posterior = analysis.posterior_sample_map['PosteriorSample-1']
        expected_paths = tuple(self.data_path(os.path.join('trees',
                'crocs-{0}.trees.gz'.format(i))) for i in range(1, 5))
        self.assertEqual(posterior.paths, expected_paths)
        self.assertEqual(posterior.burnin, 10)
        self.assertEqual(len(posterior.tip_subsets), 7)

        ts = posterior.tip_subset_map['crocodylus']
        self.assertEqual(ts.name, 'crocodylus')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 12)
        self.assertEqual(ts.tips, tuple(['poro', 'palu', 'siam', 'acut',
            'inte', 'rhom', 'more', 'nil1', 'nil2', 'john', 'nova', 'mind']))

        ts = posterior.tip_subset_map['west_niloticus']
        self.assertEqual(ts.name, 'west_niloticus')
        self.assertTrue(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 1)
        self.assertEqual(ts.tips, tuple(['nil2']))

        ts = posterior.tip_subset_map['paleosuchus']
        self.assertEqual(ts.name, 'paleosuchus')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['Ppal', 'Ptrig']))

        ts = posterior.tip_subset_map['osteolaemus']
        self.assertEqual(ts.name, 'osteolaemus')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['oste1', 'oste2']))

        ts = posterior.tip_subset_map['gharials']
        self.assertEqual(ts.name, 'gharials')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['Gav', 'Tom']))

        ts = posterior.tip_subset_map['melanosuchus']
        self.assertEqual(ts.name, 'melanosuchus')
        self.assertTrue(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 1)
        self.assertEqual(ts.tips, tuple(['Mnig']))

        ts = posterior.tip_subset_map['alligator']
        self.assertEqual(ts.name, 'alligator')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['Amis', 'Asin']))

        posterior = analysis.posterior_sample_map['PosteriorSample-2']
        expected_paths = tuple(self.data_path(os.path.join('trees',
                'gekko-{0}.trees.gz'.format(i))) for i in range(1, 5))
        self.assertEqual(posterior.paths, expected_paths)
        self.assertEqual(posterior.burnin, 10)
        self.assertEqual(len(posterior.tip_subsets), 4)

        ts = posterior.tip_subset_map['mindorensis']
        self.assertEqual(ts.name, 'mindorensis')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 15)
        self.assertEqual(ts.tips, tuple(
                'mi{0}'.format(i) for i in range(1, 16)))

        ts = posterior.tip_subset_map['negros-panay']
        self.assertEqual(ts.name, 'negros-panay')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['mi8', 'mi9']))

        ts = posterior.tip_subset_map['mindoro-caluya']
        self.assertEqual(ts.name, 'mindoro-caluya')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['mi14', 'mi15']))

        ts = posterior.tip_subset_map['kikuchii']
        self.assertEqual(ts.name, 'kikuchii')
        self.assertTrue(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 1)
        self.assertEqual(ts.tips, tuple(['mi3']))

    def test_node_age_extraction(self):
        expected_gharials = sorted([
                24.110664117758816, 23.8698119130743,
                17.525196673650544, 17.278507106852423,
                29.289878591287646, 28.736445592848177,
                22.73014866038517, 22.862970224653782,
                ])
        expected_melanosuchus = sorted([
                13.711715908355902, 13.609008882171471,
                11.853850104040147, 12.009572143207354,
                12.554953862782313, 12.812054853892633,
                14.721849319453726, 14.62532687808455,
                ])
        expected_paleosuchus = sorted([
                11.754200065064312, 11.872802801034062,
                8.50701081812682, 8.37934540330511,
                9.96578859714118, 10.157718954265865,
                9.592638399609369, 9.524296161442503,
                ])
        expected_west_niloticus = sorted([
                7.386885895079755, 7.534701038569063,
                6.269231032209964, 6.351588867999811,
                8.880262231537047, 9.049091134769006,
                5.980145075454821, 5.91632654653123,
                ])
        expected_alligator = sorted([
                57.84415042497245, 57.398644613783794,
                54.67453512303167, 55.22621447318503,
                38.76479499095923, 39.274041596147335,
                39.536140588112424, 40.32631221208484,
                ])
        expected_crocodylus = sorted([
                13.82994415, 13.72342828,
                12.243275177, 12.479067473,
                8.57111211, 8.665404431,
                11.224377229, 11.337633924
                ])
        expected_osteolaemus = sorted([
                5.094794597113713, 5.2060168284681225,
                6.6974853173279705, 6.610642417475237,
                7.621088625461037, 7.715969598889821,
                7.223120239867402, 7.0263267708481925,
                ])
        expected_kikuchii = sorted([
                0.02721530343176, 0.0261044431547,
                0.041167646421379996, 0.04907320501923,
                0.146814958946, 0.11886249552520001,
                0.1393305421926, 0.12840039221740002,
                ])
        expected_negros_panay = sorted([
                1.265195722159, 1.4052596825099999,
                0.80855471807, 0.88277607765,
                1.6664470532750002, 1.8156460958379999,
                1.462959885987, 1.491361934884,
                ])
        expected_mindoro_caluya = sorted([
                0.3925895072452, 0.3258940306721,
                0.2813710152072, 0.6968195707754,
                0.43292023300669996, 0.6484198072188,
                0.4317073941498, 0.4930305649664,
                ])
        expected_mindorensis = sorted([
                5.42999295194, 5.59535942059,
                10.11978724353, 8.99089000045,
                10.55270053408, 11.16206126123,
                6.51883087116, 7.25195430298,
                ])

        config_path = self.data_path('config-burnin98.yml')

        analysis = AnalysisManager(
                config_path = config_path,
                num_processors = 8)
        self.assertEqual(analysis.num_processors, 8)
        self.assertEqual(len(analysis.posterior_samples), 2)
        self.assertEqual(analysis.shared_nodes, [])

        analysis._extract_node_ages()

        self.assertEqual(analysis.shared_nodes, [])

        posterior = analysis.posterior_sample_map['PosteriorSample-1']
        expected_paths = tuple(self.data_path(os.path.join('trees',
                'crocs-{0}.trees.gz'.format(i))) for i in range(1, 5))
        self.assertEqual(posterior.paths, expected_paths)
        self.assertEqual(posterior.burnin, 98)
        self.assertEqual(len(posterior.tip_subsets), 7)

        ts = posterior.tip_subset_map['crocodylus']
        self.assertEqual(ts.name, 'crocodylus')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 12)
        self.assertEqual(ts.tips, tuple(['poro', 'palu', 'siam', 'acut',
            'inte', 'rhom', 'more', 'nil1', 'nil2', 'john', 'nova', 'mind']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_crocodylus):
            self.assertAlmostEqual(ages[i], expected, places = 7)

        ts = posterior.tip_subset_map['west_niloticus']
        self.assertEqual(ts.name, 'west_niloticus')
        self.assertTrue(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 1)
        self.assertEqual(ts.tips, tuple(['nil2']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_west_niloticus):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        ts = posterior.tip_subset_map['paleosuchus']
        self.assertEqual(ts.name, 'paleosuchus')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['Ppal', 'Ptrig']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_paleosuchus):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        ts = posterior.tip_subset_map['osteolaemus']
        self.assertEqual(ts.name, 'osteolaemus')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['oste1', 'oste2']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_osteolaemus):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        ts = posterior.tip_subset_map['gharials']
        self.assertEqual(ts.name, 'gharials')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['Gav', 'Tom']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_gharials):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        ts = posterior.tip_subset_map['melanosuchus']
        self.assertEqual(ts.name, 'melanosuchus')
        self.assertTrue(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 1)
        self.assertEqual(ts.tips, tuple(['Mnig']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_melanosuchus):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        ts = posterior.tip_subset_map['alligator']
        self.assertEqual(ts.name, 'alligator')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['Amis', 'Asin']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_alligator):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        posterior = analysis.posterior_sample_map['PosteriorSample-2']
        expected_paths = tuple(self.data_path(os.path.join('trees',
                'gekko-{0}.trees.gz'.format(i))) for i in range(1, 5))
        self.assertEqual(posterior.paths, expected_paths)
        self.assertEqual(posterior.burnin, 98)
        self.assertEqual(len(posterior.tip_subsets), 4)

        ts = posterior.tip_subset_map['mindorensis']
        self.assertEqual(ts.name, 'mindorensis')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 15)
        self.assertEqual(ts.tips, tuple(
                'mi{0}'.format(i) for i in range(1, 16)))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_mindorensis):
            self.assertAlmostEqual(ages[i], expected, places = 8)

        ts = posterior.tip_subset_map['negros-panay']
        self.assertEqual(ts.name, 'negros-panay')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['mi8', 'mi9']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_negros_panay):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        ts = posterior.tip_subset_map['mindoro-caluya']
        self.assertEqual(ts.name, 'mindoro-caluya')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 2)
        self.assertEqual(ts.tips, tuple(['mi14', 'mi15']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_mindoro_caluya):
            self.assertAlmostEqual(ages[i], expected, places = 10)

        ts = posterior.tip_subset_map['kikuchii']
        self.assertEqual(ts.name, 'kikuchii')
        self.assertTrue(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 1)
        self.assertEqual(ts.tips, tuple(['mi3']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_kikuchii):
            self.assertAlmostEqual(ages[i], expected, places = 10)

    def test_shared_nodes(self):
        expected_crocodylus = sorted([
                13.82994415, 13.72342828,
                12.243275177, 12.479067473,
                8.57111211, 8.665404431,
                11.224377229, 11.337633924
                ])

        config_path = self.data_path('config-burnin98-shared-nodes.yml')

        analysis = AnalysisManager(
                config_path = config_path,
                num_processors = 8)
        self.assertEqual(analysis.num_processors, 8)
        self.assertEqual(len(analysis.posterior_samples), 2)
        self.assertEqual(analysis.shared_nodes, [])

        # analysis._extract_node_ages()
        analysis.run_analysis()

        self.assertEqual(len(analysis.shared_nodes), 8)

        posterior = analysis.posterior_sample_map['PosteriorSample-1']
        expected_paths = tuple(self.data_path(os.path.join('trees',
                'crocs-{0}.trees.gz'.format(i))) for i in range(1, 5))
        self.assertEqual(posterior.paths, expected_paths)
        self.assertEqual(posterior.burnin, 98)
        self.assertEqual(len(posterior.tip_subsets), 8)

        ts = posterior.tip_subset_map['crocodylus']
        self.assertEqual(ts.name, 'crocodylus')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 12)
        self.assertEqual(ts.tips, tuple(['poro', 'palu', 'siam', 'acut',
            'inte', 'rhom', 'more', 'nil1', 'nil2', 'john', 'nova', 'mind']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_crocodylus):
            self.assertAlmostEqual(ages[i], expected, places = 7)

        ts = posterior.tip_subset_map['crocs']
        self.assertEqual(ts.name, 'crocs')
        self.assertFalse(ts.stem_based)
        self.assertIsInstance(ts.stem_based, bool)
        self.assertEqual(len(ts.tips), 12)
        self.assertEqual(ts.tips, tuple(['poro', 'palu', 'siam', 'acut',
            'inte', 'rhom', 'more', 'nil1', 'nil2', 'john', 'nova', 'mind']))
        self.assertEqual(len(ts.node_ages), 8)
        ages = sorted(ts.node_ages)
        for i, expected in enumerate(expected_crocodylus):
            self.assertAlmostEqual(ages[i], expected, places = 7)

        for shared_node in analysis.shared_nodes:
            self.assertIsInstance(shared_node.tip_subset_names, tuple)
            self.assertEqual(len(shared_node.tip_subset_names), 2)
            self.assertEqual(shared_node.tip_subset_names,
                    tuple(['crocodylus', 'crocs']))
            self.assertTrue(shared_node.tree_index in [98, 99])
            self.assertTrue(shared_node.tree_label in ['STATE_199900000',
                    'STATE_200000000'])
            path_pattern = re.compile(r'crocs-[1234].trees.gz')
            self.assertTrue(path_pattern.match(os.path.basename(
                    shared_node.tree_path)))


if __name__ == "__main__":
    if "--run-tests" in sys.argv:
        import doctest

        # doctest.testmod(verbose = True)
        suite = unittest.TestSuite()
        suite.addTest(doctest.DocTestSuite())

        tests = unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(__file__)[0])
        suite.addTests(tests)

        runner = unittest.TextTestRunner(verbosity = 2)
        runner.run(suite)

        sys.exit(0)

    main_cli()
