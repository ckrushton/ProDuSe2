import ruamel.yaml
import itertools
def generate(functions: list):
    params = ruamel.yaml.YAMLObject
    for func in functions:
        params[func.__name__]

def ValidateArgs(*validators: [function], **kwvalidators: {str:function}):
    def wrap(f):
        mergedValidators = dict(zip(f.__code__.co_varnames[:f.__code__.co_argcount], validators))
        mergedValidators.update(kwvalidators)
        def wrapped_f(*args, **kwargs):
            mergedArgs = dict(zip(f.__code__.co_varnames[:f.__code__.co_argcount], args))
            mergedArgs.update(kwargs)
            for k, v in mergedValidators:
                valArgs = {}
                if 'func' in validators[k].__code__.co_varnames[:f.__code__.co_argcount]:
                    valArgs["func"] = f
                if 'func_args' in validators[k].__code__.co_varnames[:f.__code__.co_argcount]:
                    valArgs["func_args"] = mergedArgs
                if validators[k](v, **valArgs) == False:
                    raise ValueError("{} was passed an invalid argument: {} = {}".format(f.__name__, k, v))
            f(*args, **kwargs)
        wrapped_f.validators = mergedValidators
        return wrapped_f
    return wrap

def Assert(condition: function, exception: BaseException):
    def validator(val, func, func_args):
        if not condition(val):
            if len(exception.args) and isinstance(exception.args[0], str):
                exception.args[0].format(func=func, val=val, args=func_args)
            raise exception
    return validator