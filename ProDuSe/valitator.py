import types
import os
import typing
from typing import NewType
from inspect import Signature

# If calling directly
try:
    from tator import normaliseArgs
except ImportError:
    from ProDuSe.tator import normaliseArgs

def validate(f, args: dict):
    for param, param_type in f.__annotations__.items():
        if (hasattr(param_type, '__validator__')):
            validator = param_type.__validator__
            valArgs = {}
            if 'func' in validator.__code__.co_varnames[:f.__code__.co_argcount]:
                valArgs["func"] = f
            if 'func_args' in validator.__code__.co_varnames[:f.__code__.co_argcount]:
                valArgs["func_args"] = args
            if validator(args.get(param, None), **valArgs) == False:
                if not param in args:
                    raise ValueError("Missing required argument for {}: {}".format(f.__name__, param))
                else:
                    raise ValueError("{} was passed an invalid argument: {} = {}".format(f.__name__, param, args.get(param, None)))

def Validate():
    def wrap(f):
        def wrapped_f(*args, **kwargs):
            validate(f, normaliseArgs(args, kwargs), wrapped_f.validators)
            f(*args, **kwargs)
        return wrapped_f
    return wrap

def Assert(condition: types.FunctionType, exception: BaseException):
    def validator(val, func, func_args):
        if not condition(val):
            if len(exception.args) and isinstance(exception.args[0], str):
                exception.args[0].format(func=func, val=val, args=func_args)
            raise exception
    return validator

_T = typing.TypeVar('_T')
def ValidatedType(name: str, tp: typing.Type[_T], validator) -> typing.Type[_T]:
    """
    Extends the functionality of typing.NewType to attach the validation function
    :param name: See typing.NewType(name)
    :param tp:   See typing.NewType(name)
    :param validator: The validation function called to validate the value of the variable
    :return: See typing.NewType
    """
    x = NewType(name, tp)
    x.__validator__ = validator
    return x

PathOrNone = ValidatedType('PathOrNone', str, lambda x: not x or os.path.isfile(x))
Path = ValidatedType('Path', str, os.path.isfile)
UnsignedInt = ValidatedType('UnsignedInt', int, lambda x: x >= 0)

def Domain(min, max) -> type:
    """
    Create a Type restricted to the specified domain
    :param min: The lower bound of the domain
    :param max: The upper bound of the domain
    :return: A ValidatedType representing the specified domain
    """
    return ValidatedType('Domain_{}-{}'.format(min, max), int, lambda x: min <= x <= max)

def which(program: str):
    """
    Locates the true path to the executable in the environment
    Borrowed from https://stackoverflow.com/a/377028
    :param program: The program to determine a path to
    :return: The fully qualified path to the executable or None if not found
    """
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    tmp = fname.split(' ')
    if len(tmp) > 1: #Handle path containing arguments
        fname = tmp[0]
    if fpath:
        if is_exe(os.path.join(fpath, fname)):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, fpath, fname)
            if is_exe(exe_file):
                return exe_file

    return None

Executable = ValidatedType('Executable', str, which)
