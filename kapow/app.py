# Copyright (C) 2018 Christopher Cave-Ayland
# You may use, distribute and modify this code under the terms of the
# MIT open source license. See https://opensource.org/licenses/MIT if
# a copy was not provided with this code.

from collections import defaultdict
from functools import update_wrapper
import inspect
import sys
from simtk.openmm import app
from .class_map import class_map
from .interface import Interface
from .wrappers import wrap_method

exclusions = defaultdict(lambda: {}, [])

skip = ['CharmmPSFWarning']

module = sys.modules[__name__]
for name in dir(app):
    if name.startswith('_') or name in skip:
        continue
    attr = getattr(app, name)
    if inspect.isclass(attr):
        # does this logic cause problems for reloads?
        if getattr(module, name, None) is None:
            openmm_class = attr
            new_class = Interface(
                name,
                [openmm_class],
                {'exclude': exclusions[name], 'preserve': []})
            setattr(module, name, new_class)
        class_map[attr] = getattr(module, name)
    elif inspect.isfunction(attr):
        # doesn't seem to actually get used
        new_func = wrap_method(attr)
        update_wrapper(new_func, attr)
        setattr(module, name, update_wrapper(wrap_method, attr))
    else:
        setattr(module, name, attr)
