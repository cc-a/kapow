import inspect
from .class_map import class_map
from .help_strings import class_help_string
from .wrappers import build_ArrayWrapper, build_DictWrapper
from .wrappers import ValueWrapper, wrap_method


def create_init_function(cls):
    @wrap_method
    def __init__(self, *args, **kwargs):
        if kwargs.get('create_wrapped', True):
            self.wrapped_object = cls(*args, **kwargs)
        else:
            self.wrapped_object = None
    return __init__


def create_wrap_function(cls):
    def Wrap(inst):
        wrapper = class_map[cls](create_wrapped=False)
        wrapper.wrapped_object = inst
        return wrapper
    return Wrap


class Pythonize(type):
    """Metaclass to provide a more Pythonic interface to openmm objects.

    Where values are controlled through getter and setter routines a single
    property interface is created. Where getters, setters and optionally
    add and remove routines are used to interact with array elements a
    more sophisticated interface is created providing List like semantics.
    """

    def __new__(mcs, name, bases, attrs):

        if '__init__' not in attrs:
            attrs['__init__'] = create_init_function(bases[0])
        attrs['Wrap'] = staticmethod(create_wrap_function(bases[0]))

        attrs_to_add = {}
        for base_class in bases[0].__mro__[:-1]:
            attrs_to_add.update(base_class.__dict__)

        for method in attrs['exclude']:
            attrs_to_add.pop(method)
        for method in attrs['preserve']:
            attrs_to_add.pop(method)
            attrs[method] = getattr(bases[0], method)

        all_array_roots = [attr[6:-1]
                           for attr in attrs_to_add.keys()
                           if 'getNum' in attr]

        array_getters = {}
        for root in all_array_roots:
            array_getters[root] = [
                attr for attr in attrs_to_add.keys()
                if 'get' + root in attr]

        array_setters = {}
        for root in all_array_roots:
            array_setters[root] = [attr for attr in attrs_to_add.keys()
                                   if 'set' + root in attr]

        array_adders = {}
        for root in all_array_roots:
            try:
                array_adders[root] = [attr for attr in attrs_to_add.keys()
                                      if 'add' + root in attr][0]
            except IndexError:
                array_adders[root] = None

        array_removers = {}
        for root in all_array_roots:
            try:
                array_removers[root] = [attr for attr in attrs_to_add.keys()
                                        if 'remove' + root in attr][0]
            except IndexError:
                array_removers[root] = None

        for root in all_array_roots:
            method_name = root[0].lower() + root[1:] + 's'
            try:
                args = ('getNum%ss' % root, array_getters[root],
                        array_setters[root], array_adders[root],
                        array_removers[root])
                attrs[method_name] = build_ArrayWrapper(method_name,
                                                        bases[0], *args)
            except ValueError:
                # no getters were found that go along with this array
                pass

        used_methods = set(method for method_type in (
                                                      array_getters,
                                                      array_setters)
                           for methods in method_type.values()
                           for method in methods)
        used_methods.update(array_adders.values())
        used_methods.update(array_removers.values())
        used_methods.update('getNum%ss' % root for root in all_array_roots)

        for attr in attrs_to_add:
            if attr.startswith('get') and attr not in used_methods:
                if attr[4:5].islower():
                    method_name = attr[3:4].lower() + attr[4:]
                else:
                    method_name = attr[3:]

                setter = attr.replace('get', 'set')
                if setter in attrs_to_add:
                    try:
                        attrs[method_name] = ValueWrapper(
                            bases[0],
                            attr,
                            setter)
                        used_methods.update((attr, setter))
                    except TypeError:
                        arg_spec = inspect.getargspec(
                            getattr(bases[0], attr))
                        if arg_spec.defaults is None or \
                           len(arg_spec.args) == 2:
                            dw = build_DictWrapper(
                                method_name, getattr(bases[0], attr),
                                getattr(bases[0], setter))
                            attrs[method_name + 's'] = dw
                else:
                    try:
                        attrs[method_name] = ValueWrapper(
                            bases[0],
                            attr)
                        used_methods.add(attr)
                    except TypeError:
                        # deals with some special cases such as
                        # NonbondedForce.getPMEParametersInContext
                        pass

        for attr in attrs_to_add:
            if attr not in used_methods and not attr.startswith('__'):
                attrs[attr] = attrs_to_add[attr]

        attrs['__doc__'] = class_help_string.format(bases[0], bases[0].__doc__)
        attrs.pop('exclude')
        attrs.pop('preserve')
        return type.__new__(mcs, name, (), attrs)
