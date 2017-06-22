import re

from tator import normaliseArgs

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap
import jmespath

argDocRegEx = re.compile(r':param ([^:]+):\s*(.+)$', flags=re.MULTILINE)

def generate(functions: list):
    config = {}
    for func in functions:
        params = CommentedMap()
        params.update({arg:None for arg in func.__code__.co_varnames[:func.__code__.co_argcount]})
        for match in argDocRegEx.finditer(func.__doc__):
            params.yaml_add_eol_comment(match.group(2), match.group(1))
        config[func.__name__] = params
    return config

def map(func, yaml_node):
    args = {}
    for arg, expression in func.configMap:
        val = expression.search(yaml_node)
        if val is None:
            val = yaml_node.get(arg, None)
        args[arg] = val
    return args

def ConfigMap(*args, **kwargs):
    def wrap(f):
        f.configMap = {arg:jmespath.compile(arg) for arg in f.__code__.co_varnames[:f.__code__.co_argcount]}
        for arg, expression in normaliseArgs(args, kwargs):
            f.configMap[arg] = jmespath.compile(expression)
        return f
    return wrap

def multidoc_anchor_composer(self):
    """
    https://stackoverflow.com/a/40702117
    :param self:
    :return: The composed node
    """
    self.get_event()
    node = self.compose_node(None, None)
    self.get_event()
    return node

def construct_concat(loader, node):
        return "".join([str(i) for i in loader.construct_sequence(node)])

ruamel.yaml.composer.Composer.compose_document = multidoc_anchor_composer
ruamel.yaml.add_constructor(u'!concat', construct_concat)