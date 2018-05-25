# Copyright (C) 2018 Christopher Cave-Ayland
# You may use, distribute and modify this code under the terms of the
# MIT open source license. See https://opensource.org/licenses/MIT if
# a copy was not provided with this code.

from functools import update_wrapper
import inspect
from .wrappers import wrap_method, ClassMethodWrapper


class Interface(type):
    """Metaclass to provide a more Pythonic interface to openmm objects.

    Where values are controlled through getter and setter routines a single
    property interface is created. Where getters, setters and optionally
    add and remove routines are used to interact with array elements a
    more sophisticated interface is created providing List like semantics.
    """

    def __new__(mcs, name, bases, attrs):
        attrs_to_add = {}
        for base_class in bases[0].__mro__[:-1]:
            attrs_to_add.update(base_class.__dict__)

        for attr_name in attrs_to_add:
            attr = attrs_to_add[attr_name]
            if attr_name.startswith('_'):
                # special methods pass unchanged
                attrs[attr_name] = attr
            elif inspect.isfunction(attr):
                # functions are wrapped to handle adapter classes
                wrapped_attr = wrap_method(attr)
                update_wrapper(wrapped_attr, attr)
                attrs[attr_name] = wrapped_attr
            elif isinstance(attr, classmethod):
                # classmethods need special treatment
                attrs[attr_name] = ClassMethodWrapper(attr)
            else:
                # anything else just pass it through
                # a lot of staticmethods go through here
                # that probably need wrapping?
                attrs[attr_name] = attr

        new_type = type.__new__(mcs, name, (), attrs)
        return new_type
