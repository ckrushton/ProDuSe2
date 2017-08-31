import types

def normaliseArgs(func: types.FunctionType, args: list, kwargs: dict) -> dict:
    mergedArgs = dict(zip(func.__code__.co_varnames[:func.__code__.co_argcount], args))
    mergedArgs.update(kwargs)
    return mergedArgs