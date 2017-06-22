import types
import os
from typing import TypeVar

from tator import normaliseArgs

def validate(f, args):
    if not hasattr(f, 'validators'):
        raise RuntimeError("The passed function needs to have a @ValidateArgs annotation to use this function.")
    for k, v in f.validators:
        valArgs = {}
        if 'func' in f.validators[k].__code__.co_varnames[:f.__code__.co_argcount]:
            valArgs["func"] = f
        if 'func_args' in f.validators[k].__code__.co_varnames[:f.__code__.co_argcount]:
            valArgs["func_args"] = args
        if f.validators[k](v, **valArgs) == False:
            raise ValueError("{} was passed an invalid argument: {} = {}".format(f.__name__, k, v))

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

class Type:
    def __new__(cls, t: type, validator: types.FunctionType):
        x = t.__new__(t)
        x.__validator__ = validator
        return x

PathOrNone = Type(str, lambda x: not x or os.path.isfile(x))
Path = Type(str, os.path.isfile)

UnsignedInt = Type(int, lambda x: x >= 0)

def Domain(min, max):
    return Type(int, lambda x: min <= x <= max)