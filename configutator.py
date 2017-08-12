import re
import os
from getopt import gnu_getopt
from asciimatics.screen import Screen, StopApplication, NextScene
from asciimatics.scene import Scene
import asciimatics.widgets as Widgets
from inspect import signature, getdoc, Parameter
from distutils.util import strtobool

from tator import normaliseArgs
from valitator import validate

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap
import jmespath

YAMLLoader = ruamel.yaml.SafeLoader
Widgets.Frame.palette['edit_text'] = (7,2,7)

funcDocRegex = re.compile(r'^[^:]+?(?=[\s]*:[\w]+)')
argDocRegEx = re.compile(r':param ([^:]+):[\t ]*(.+)$', flags=re.MULTILINE)
#cmdlineRegEx = re.compile(r' *(?:-+([^= \'\"]+)[= ]?)?(?:([\'\"])([^\2]+?)\2|([^- \"\']+))?')
#
#def parseArgs(argv):
#    return [(match.group(1), match.group(3)+match.group(4)) for match in cmdlineRegEx.finditer("".join(argv))]

def getParamDoc(func) -> dict:
    return {match.group(1): match.group(2) for match in argDocRegEx.finditer(getdoc(func) or "")}

def generate(functions: list):
    config = {}
    for func in functions:
        sig = signature(func)
        params = CommentedMap()
        params.update({param:None for param, _ in sig.parameters})
        for match in argDocRegEx.finditer(getdoc(func)):
            params.yaml_add_eol_comment(match.group(2), match.group(1))
        config[func.__name__] = params
    return config

def map(func, yaml_node, path = None) -> dict:
    if path:
        yaml_node = yaml_node.get(path, yaml_node)
    elif hasattr(func, '__cfgMap__') and '_func' in func.__cfgMap__:
        yaml_node = func.__cfgMap__['_func'].search(yaml_node)
    args = {}
    sig = signature(func)
    for name, param in sig.parameters.items(): #type: str, Parameter
        val = yaml_node.get(name, param.default)
        if val != param.empty:
            args[name] = val
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
        sig = signature(f)
        if not hasattr(f, '__cfgMap__'):
            f.__cfgMap__ = {param:jmespath.compile(param) for param in sig.parameters}
        for arg, expression in normaliseArgs(f, args, kwargs).items():
            if arg == '_func' or not expression:
                f.__cfgMap__[arg] = expression
            else:
                f.__cfgMap__[arg] = jmespath.compile(expression)
        return f
    return wrap

def configUnmapped(func, param: str):
    return hasattr(func, '__cfgMap__') and param in func.__cfgMap__ and func.__cfgMap__[param] is None

def boolParam(param: Parameter, desc: str, layout: Widgets.Layout):
    layout.add_widget(Widgets.CheckBox(desc, param.name, param.name))

def textParam(param: Parameter, desc: str, layout: Widgets.Layout):
    layout.add_widget(Widgets.Text(name=param.name, label=param.name))
    if desc:
        layout.add_widget(Widgets.Label(desc))
        layout.add_widget(Widgets.Divider(False))

widgetMap = {
    'Path': textParam,
    'PathOrNone': textParam,
    str: textParam,
    int: textParam,
    bool: boolParam
}

def configUI(screen, functions: list, yaml_node, title=''):
    def _go():
        raise StopApplication('')
    def _save():
        pass
    def _load():
        pass
    window = Widgets.Frame(screen, screen.height * 2 // 3, screen.width * 2 // 3, title=title)
    scene = Scene([window])

    for f in functions:
        doc = getParamDoc(f)
        layout = Widgets.Layout([1], True)
        window.add_layout(layout)
        sig = signature(f)
        for name, param in sig.parameters.items(): #type: str, Parameter
            if param.annotation in widgetMap:
                widgetMap[param.annotation](param, doc.get(name, ''), layout)

    toolbar = Widgets.Layout([1,1,1])
    window.add_layout(toolbar)
    load_button = Widgets.Button('Load', _load)
    save_button = Widgets.Button('Save', _save)
    go_button = Widgets.Button('GO', _go)
    toolbar.add_widget(load_button, 0)
    toolbar.add_widget(save_button, 1)
    toolbar.add_widget(go_button, 2)

    window.fix()
    screen.play([scene])

def mapArgs(optList, functions) -> (dict, dict):
    overrides = {}
    for f in functions:
        sig = signature(f)
        overrides[f] = {}
        for opt, val in optList:
            if opt == '--'+f.__name__+':':
                ConfigMap(_func=val)(f)
            else:
                name = opt[len(f.__name__) + 3 if opt.startswith('--' + f.__name__ + ':') else 2:]
                if name == 'config' or name.startswith('config:'): continue
                if sig.parameters[name].annotation is bool: val = strtobool(val or 'True')
                elif sig.parameters[name].annotation is not Parameter.empty: val = sig.parameters[name].annotation(val)
                if opt[-1] == ':':
                    ConfigMap(**{name:val})(f)
                else:
                    overrides[f][name] = val
    return overrides

def loadConfig(args: list, functions: tuple, configParam='config', defaultConfig=None, title='') -> (dict, list):
    call_name = args[0]
    args = args[1:]
    options = [configParam + '=', '{}:params:='.format(configParam), '{}:batch:='.format(configParam)]
    for f in functions:
        sig = signature(f)
        if len(functions) == 1:
            options += [name + ('' if param.annotation is bool and param.default == False else '=') for name, param in sig.parameters.items() if not configUnmapped(f, name)]
        else:
            options += [f.__name__ + ':' + param + '=' for param in sig.parameters if not configUnmapped(f, param)]
    optList, originalParams = gnu_getopt(args, 'h', options)
    params = originalParams.copy()

    # Search for help switch in passed parameters
    for opt, val in optList:
        if opt == '-h':
            from sys import stderr
            if title: stderr.write(title + '\n')
            for f in functions:
                sig = signature(f)
                if len(functions) > 1: stderr.write("{}:\n".format(f.__name__))
                funcDocs = funcDocRegex.findall(getdoc(f))
                if len(funcDocs):
                    stderr.write(funcDocs[0] + '\n')
                doc = getParamDoc(f)
                for name, param in sig.parameters.items(): #type:str, Parameter
                    if hasattr(f, '__cfgMap__') and f.__cfgMap__.get(name, None) is not None:
                        argName = (f.__name__ + ':' if len(functions) > 1 else '') + name
                        default = ''
                        if param.default != param.empty:
                            default = "[Default: {}]".format(param.default)
                        stderr.write("  --{} - {} {}\n".format(argName, doc.get(name, ''), default))
            stderr.write("  --{} - Path to configuration file. All other specified options override config. " 
                         "If no config file specified, all options without default must be specified.\n".format(configParam))
            exit()

    # Search for configParam in passed parameters
    for opt, val in optList:
        if opt == '--' + configParam:
            if not (os.path.isfile(val) and os.access(val, os.R_OK)):
                raise FileNotFoundError("The passed config path was not found: {}".format(val))
            configExpression = None
            paramExpression = None
            batchExpression = None
            for o, v in optList:
                if o == '--{}:'.format(configParam):
                    configExpression = jmespath.compile(v)
                elif o == '--{}:params:'.format(configParam):
                    paramExpression = jmespath.compile(v)
                elif o == '--{}:batch:'.format(configParam):
                    batchExpression = jmespath.compile(v)
            overrides = mapArgs(optList, functions)
            cfgs = ruamel.yaml.load_all(open(val), YAMLLoader)
            for cfg in cfgs:
                if configExpression: cfg = configExpression.search(cfg)
                if batchExpression:
                    jobs = batchExpression.search(cfg)
                else:
                    jobs = [cfg]
                for job in jobs:
                    if paramExpression:
                        params = originalParams.copy()
                        params += paramExpression.search(job)
                    argMap = {}
                    for f in functions:
                        mapping = map(f, job)
                        argMap[f] = mapping
                        validate(f, mapping)
                    for f, m in overrides.items():
                        if f not in argMap:
                            argMap[f] = {}
                        argMap[f].update(m)
                    yield argMap, params
            return

    if len(args):
        # If parameters were passed, skip TUI
        yield mapArgs(optList, functions)
        return
    else:
        # No parameters were passed, load TUI
        cfg = None
        if defaultConfig is not None:
            cfg = ruamel.yaml.load(defaultConfig, YAMLLoader)
        Screen.wrapper(configUI, arguments=(functions, cfg, title))
