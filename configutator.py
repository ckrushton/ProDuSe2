import re
import os
from getopt import gnu_getopt
from asciimatics.screen import Screen
from time import sleep

from tator import normaliseArgs
from valitator import validate

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap
import jmespath

argDocRegEx = re.compile(r':param ([^:]+):\s*(.+)$', flags=re.MULTILINE)
cmdLineRegEx = re.compile(r'')

def generate(functions: list):
    config = {}
    for func in functions:
        params = CommentedMap()
        params.update({arg:None for arg in func.__code__.co_varnames[:func.__code__.co_argcount]})
        for match in argDocRegEx.finditer(func.__doc__):
            params.yaml_add_eol_comment(match.group(2), match.group(1))
        config[func.__name__] = params
    return config

def map(func, yaml_node) -> list:
    if hasattr(func, '__cfgMap__') and '_func' in func.__cfgMap__:
        yaml_node = yaml_node.get(func.__cfgMap__['_func'], yaml_node)
    args = {}
    for arg in func.__code__.co_varnames[:func.__code__.co_argcount]:
        val = yaml_node.get(arg, None)
        if val:
            args[arg] = val
    if hasattr(func, '__cfgMap__'):
        for arg, expression in func.__cfgMap__.items():
            if arg != '_func' and expression:
                val = expression.search(yaml_node)
                if val is None:
                    val = yaml_node.get(arg, None)
                args[arg] = val
    return args

def ConfigMap(*args, **kwargs):
    def wrap(f):
        if not hasattr(f, '__cfgMap__'):
            f.__cfgMap__ = {arg:jmespath.compile(arg) for arg in f.__code__.co_varnames[:f.__code__.co_argcount]}
        for arg, expression in normaliseArgs(f, args, kwargs).items():
            if arg == '_func' or not expression:
                f.__cfgMap__[arg] = expression
            else:
                f.__cfgMap__[arg] = jmespath.compile(expression)
        return f
    return wrap

def construct_concat(loader, node):
    return "".join([str(i) for i in loader.construct_sequence(node)])

def construct_ignored(loader: ruamel.yaml.BaseLoader, node):
    obj = loader.construct_mapping(node) # type: object
    if not hasattr(obj, '__dict__'):
        obj = type(type(obj).__name__ + '_ignored', (type(obj),), {'__ignore__':True})(obj)
    else:
        obj.__ignore__ = True
    return obj

class StreamAnchorComposer():
    def compose_document(self):
        self.get_event()
        node = self.compose_node(None, None)
        self.get_event()
        return node

class BaseLoader(StreamAnchorComposer, ruamel.yaml.BaseLoader): pass
class SafeLoader(StreamAnchorComposer, ruamel.yaml.SafeLoader): pass
class Loader(StreamAnchorComposer, ruamel.yaml.Loader): pass
class RoundTripLoader(StreamAnchorComposer, ruamel.yaml.RoundTripLoader): pass

YAMLLoader = Loader

YAMLLoader.add_constructor(u'!concat', construct_concat)
YAMLLoader.add_constructor(u'!ignore', construct_ignored)

def handleArgs(args: list, functions: tuple, configParam = '--config'):
    optlist, params = gnu_getopt(args, '', ['config=',])
    for opt, val in optlist:
        if opt == '--config':
            if not (os.path.isfile(val) and os.access(val, os.R_OK)):
                raise FileNotFoundError("The passed config path was not found: {}".format(val))
            cfgs = ruamel.yaml.load_all(open(val), YAMLLoader)
            for cfg in cfgs:
                argMap = {}
                if not hasattr(cfg, '__ignore__'):
                    for f in functions:
                        mapping = map(f, cfg)
                        argMap[f] = mapping
                        validate(f, mapping)
                    yield argMap
        # else:
        #     screen = Screen.open(3)
        #     screen.print_at("Hello", 0, 0)
        #     screen.refresh()
        #     sleep(10)
        #     screen.close()

