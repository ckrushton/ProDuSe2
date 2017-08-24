import re, os, pickle, sys
from getopt import gnu_getopt
from inspect import signature, getdoc, getfile, Parameter

# If calling directly, this works fine
try:
    from tator import normaliseArgs
    from valitator import validate
# If installed
except ModuleNotFoundError:
    from ProDuSe.tator import normaliseArgs
    from ProDuSe.valitator import validate

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap
import jmespath

from asciimatics.screen import Screen, StopApplication, NextScene
from asciimatics.scene import Scene
import asciimatics.widgets as Widgets
from asciimatics.event import KeyboardEvent

YAMLLoader = ruamel.yaml.SafeLoader

funcDocRegex = re.compile(r'^[^:]+?(?=[\s]*:[\w]+)')
argDocRegEx = re.compile(r':param ([^:]+):[\t ]*(.+)$', flags=re.MULTILINE)

def getParamDoc(func) -> dict:
    return {match.group(1): match.group(2) for match in argDocRegEx.finditer(getdoc(func) or "")}

def getTrueAnnotationType(annotation):
    return getattr(annotation, '__supertype__', annotation)

def strtobool(val: str) -> bool:
    """Convert a string representation of truth to true or false.

    True values are 'y', 'yes', 't', 'true', 'on', and '1'.
    False values are 'n', 'no', 'f', 'false', 'off', and '0'.
    Raises ValueError if 'val' is anything else.

    Taken from distutils.util and modified to return correct type.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return True
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return False
    else:
        raise ValueError("invalid truth value %r" % (val,))

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

# -- Config file mapping --

def ConfigMap(*args, **kwargs):
    def wrap(func):
        sig = signature(func)
        nArgs = normaliseArgs(func, args, kwargs)
        if not hasattr(func, '__cfgMap__'):
            func.__cfgMap__ = {}
        if '_func' in nArgs:
            func.__cfgMap__['_func'] = jmespath.compile(nArgs['_func'])
        for param in sig.parameters:
            if param in nArgs:
                expression = nArgs.get(param)
                func.__cfgMap__[param] = jmespath.compile(expression) if expression else None
            elif param not in func.__cfgMap__:
                func.__cfgMap__[param] = jmespath.compile(param)

        return func
    return wrap

def mapConfig(func, yaml_node, path = None) -> dict:
    if path:
        yaml_node = yaml_node.get(path, yaml_node)
    elif hasattr(func, '__cfgMap__') and '_func' in func.__cfgMap__:
        yaml_node = func.__cfgMap__['_func'].search(yaml_node)
    args = {}
    sig = signature(func)
    for name, param in sig.parameters.items(): #type: str, Parameter
        if hasattr(func, '__cfgMap__') and name in func.__cfgMap__:
            expression = func.__cfgMap__[name]
            if expression:
                val = expression.search(yaml_node)
                if val is None:
                    val = param.default
                args[name] = val
        else:
            val = yaml_node.get(name, param.default)
            if val != param.empty:
                args[name] = val
    return args

def configUnmapped(func, param: str):
    return hasattr(func, '__cfgMap__') and param in func.__cfgMap__ and func.__cfgMap__[param] is None

# -- Command line argument mapping --

def ArgMap(*args, **kwargs):
    def wrap(func):
        sig = signature(func)
        nArgs = normaliseArgs(func, args, kwargs)
        if not hasattr(func, '__argMap__'):
            func.__argMap__ = {}
        if '_func' in nArgs:
            func.__argMap__['_func'] = nArgs['_func']
        for param in sig.parameters:
            if param in nArgs:
                func.__argMap__[param] = nArgs.get(param) or None
            #elif param not in func.__argMap__:
            #    func.__argMap__[param] = param

        return func
    return wrap

def mapArgs(optList, functions) -> (dict, dict):
    args = {}
    for f in functions:
        sig = signature(f)
        # Handle _func in argMap
        prefix = f.__name__ + ':'
        if hasattr(f, '__argMap__') and '_func' in f.__argMap__ and f.__argMap__['_func'] is not None:
            prefix = f.__argMap__['_func']
            if prefix != '':
                prefix += ':'
        prefix = '--' + prefix
        args[f] = {}
        for name, param in sig.parameters.items():
            argName = resolveArgName(f, name)
            if not argName: continue
            for opt, val in optList:
                if opt == '--' + argName:
                    if issubclass(getTrueAnnotationType(sig.parameters[name].annotation), bool):
                        val = strtobool(val or 'True')
                    elif sig.parameters[name].annotation is not Parameter.empty:
                        val = getTrueAnnotationType(sig.parameters[name].annotation)(val)
                    args[f][name] = val

        # Handle config path overrides TODO: This should be moved out of this function
        for opt, val in optList:
            if opt == prefix:
                ConfigMap(_func=jmespath.compile(val))(f)
            elif opt == prefix + name + ':':
                ConfigMap(**{name:jmespath.compile(val)})(f)
    return args

def resolveArgName(func, name):
    # Prefix parameter with name if argMap unspecified
    if hasattr(func, '__argMap__'):
        if name in func.__argMap__:
            if func.__argMap__[name]:
                return func.__argMap__[name]
            else:
                return None
        elif '_func' in func.__argMap__:
            if func.__argMap__['_func'] is not None:
                return func.__argMap__['_func'] + ':' + name
            else:
                return name
        else:
            return name
    else:
        return func.__name__ + ':' + name

def argUnmapped(func, param: str):
    return hasattr(func, '__argMap__') and param in func.__argMap__ and func.__argMap__[param] is None

# -- TUI Builder --

def boolParam(param: Parameter, desc: str, layout: Widgets.Layout, maxWidth = None):
    layout.add_widget(Widgets.CheckBox(desc, param.name, param.name))

def textParam(param: Parameter, desc: str, layout: Widgets.Layout, maxWidth = None):
    layout.add_widget(Widgets.Text(name=param.name, label=param.name))
    if desc:
        lastI = 0
        for i in range(maxWidth or len(desc), len(desc), maxWidth or len(desc)):
            layout.add_widget(Widgets.Label(desc[lastI:i]))
            lastI = i
        layout.add_widget(Widgets.Divider(False))

widgetMap = {
    'Path': textParam,
    'PathOrNone': textParam,
    'str': textParam,
    'int': textParam,
    'bool': boolParam
}

def configUI(screen, functions: list, yaml_node, call_name, title=''):
    Widgets.Frame.palette['edit_text'] = (Screen.COLOUR_BLACK, Screen.A_NORMAL, Screen.COLOUR_CYAN)
    Widgets.Frame.palette['section_header'] = (Screen.COLOUR_GREEN, Screen.A_UNDERLINE, Screen.COLOUR_BLUE)

    windowWidth = screen.width * 2 // 3
    data = {}
    historyPath = os.path.expanduser('~/.configutator_history')
    if os.path.exists(historyPath):
        with open(historyPath, 'rb') as file:
            history = pickle.load(file)
            if call_name in history:
                data = history[call_name]
    window = Widgets.Frame(screen, screen.height * 2 // 3, windowWidth, title=title, data=data)
    scene = Scene([window])

    def saveHistory():
        history = {}
        if os.path.exists(historyPath):
            with open(historyPath, 'rb') as file:
                history = pickle.load(file)
        history[call_name] = window.data
        with open(historyPath, 'wb+') as file:
            pickle.dump(history, file)

    def _go():
        window.save()
        saveHistory()
        raise StopApplication('')
    def _close():
        window.save()
        saveHistory()
        raise StopApplication('')
    def _save():
        pass
    def _load():
        pass
    def inputHandler(event):
        if isinstance(event, KeyboardEvent):
            c = event.key_code
            # Stop on ctrl+q or ctrl+x or esc
            if c in (17, 24, 27):
                _close()

    for f in functions:
        doc = getParamDoc(f)
        layout = Widgets.Layout([1])
        window.add_layout(layout)
        sig = signature(f)
        if len(functions) > 1:
            label = Widgets.Label(f.__name__+':')
            label.custom_colour = 'section_header'
            layout.add_widget(label)
        for name, param in sig.parameters.items(): #type: str, Parameter
            if param.annotation.__name__ in widgetMap:
                widgetMap[param.annotation.__name__](param, doc.get(name, ''), layout, windowWidth)
        layout.add_widget(Widgets.Divider(False))

    toolbar = Widgets.Layout([1,1,1,1])
    window.add_layout(toolbar)
    load_button = Widgets.Button('Load', _load)
    save_button = Widgets.Button('Save', _save)
    cancel_button = Widgets.Button('Cancel', _close)
    go_button = Widgets.Button('GO', _go)
    toolbar.add_widget(load_button, 0)
    toolbar.add_widget(save_button, 1)
    toolbar.add_widget(cancel_button, 2)
    toolbar.add_widget(go_button, 3)

    window.fix()
    screen.play([scene], unhandled_input=inputHandler)

# -- Loader --

def loadConfig(args: list, functions: tuple, title='', positionalDoc: list=[], configParam='config', defaultConfig=None, configExpression = None, paramExpression = None, batchExpression=None) -> (dict, list):
    if isinstance(configExpression, str):
        configExpression = jmespath.compile(configExpression)
    if isinstance(paramExpression, str):
        paramExpression = jmespath.compile(paramExpression)
    if isinstance(batchExpression, str):
        batchExpression = jmespath.compile(batchExpression)

    call_name = args[0]
    args = args[1:]
    options = [configParam + '=', '{}:='.format(configParam), '{}:params:='.format(configParam), '{}:batch:='.format(configParam)]
    for f in functions:
        # Build command line argument list from function parameters
        sig = signature(f)
        for name, param in sig.parameters.items():
            name = resolveArgName(f, name)
            if not name: continue
            # If boolean, make parameter valueless
            if not ((issubclass(getTrueAnnotationType(param.annotation), bool) or isinstance(param.default, bool)) and param.default == False):
                name += '='
            options.append(name)
    optList, originalParams = gnu_getopt(args, 'h', options)
    params = originalParams.copy()

    # Search for help switch in passed parameters
    for opt, val in optList:
        if opt == '-h':
            from sys import stderr
            import __main__
            try:
                if title: stderr.write("{}\tUse: {} [options] {}\n{}\n".format(
                                        title,
                                        getfile(__main__),
                                        ' '.join([name for name, _ in positionalDoc]),
                                        '\n'.join([name + ' - ' + desc for name, desc in positionalDoc if desc])
                                        ))
            except TypeError:
                pass
            for f in functions:
                # Build command line argument option descriptions from function parameters
                sig = signature(f)
                if len(functions) > 1: stderr.write("{}:\n".format(f.__name__))
                funcDocs = funcDocRegex.findall(getdoc(f))
                if len(funcDocs):
                    stderr.write(funcDocs[0] + '\n')
                doc = getParamDoc(f)
                for name, param in sig.parameters.items(): #type:str, Parameter
                    argName = resolveArgName(f, name)
                    if not argName: continue
                    default = ''
                    if param.default != param.empty:
                        default = "[Default: {}]".format(param.default)
                    stderr.write("  --{} - {} {}\n".format(argName, doc.get(name, ''), default))

            stderr.write("General:\n  --{} - Path to configuration file. All other specified options override config. " 
                         "If no config file specified, all options without default must be specified.\n".format(configParam))
            exit()

    # Search for configParam in passed parameters
    for opt, val in optList:
        if opt == '--' + configParam:
            if not (os.path.isfile(val) and os.access(val, os.R_OK)):
                raise FileNotFoundError("The config path was not found: {}".format(val))
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
                if configExpression:
                    tmp = configExpression.search(cfg)
                    if tmp != None:
                        cfg = tmp
                if batchExpression:
                    jobs = batchExpression.search(cfg)
                    if jobs == None:
                        jobs = [cfg]
                else:
                    jobs = [cfg]
                for job in jobs:
                    if paramExpression:
                        params = originalParams.copy()
                        tmp = paramExpression.search(job)
                        if isinstance(tmp, list):
                            params += tmp
                    argMap = {}
                    for f in functions:
                        mapping = mapConfig(f, job)
                        argMap[f] = mapping
                        validate(f, mapping)
                    for f, m in overrides.items():
                        if f not in argMap:
                            argMap[f] = {}
                        argMap[f].update(m)
                    #TODO os.fork() to parallelize jobs
                    yield argMap, params
            return

    if len(args):
        # If parameters were passed, skip TUI
        yield mapArgs(optList, functions), params
        return
    else:
        # No parameters were passed, load TUI
        cfg = None
        if defaultConfig is not None:
            cfg = ruamel.yaml.load(defaultConfig, YAMLLoader)
        Screen.wrapper(configUI, arguments=(functions, cfg, call_name, title))
