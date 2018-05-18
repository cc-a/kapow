import inspect
from .wrappers import ArrayWrapper, DictWrapper, ValueWrapper, build_ArrayWrapper
from class_map import class_map


def create_init_function(cls):
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
                if 'get' + root in attr]  # and
                # isinstance(bases[0].__dict__[attr], staticmethod)]
                # inspect.getargspec(getattr(bases[0], attr)).args[0] == 'self']

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
                # new_class = type(
                #     'tmpName', (ArrayWrapper,),
                #     {'__doc__': ArrayWrapper.__doc__ % args})
                attrs[method_name] = build_ArrayWrapper(method_name,
                                                        bases[0], *args)
                # attrs[method_name] = ArrayWrapper(
                #     bases[0],
                #     len_='getNum%ss' % root,
                #     getters=array_getters[root],
                #     setters=array_setters[root],
                #     adder=array_adders[root],
                #     remover=array_removers[root])
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

        # used_methods = set(method for methods in array_adders.values()
        #                    for method in methods)
        # used_methods.update(method for methods in array_removers.values()
        #                     for method in methods)
        # used_methods.update(method for methods in array_getters.values()
        #                     for method in methods)
        # used_methods.update(method for methods in array_setters.values()
        #                     for method in methods)

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
                            attrs[method_name + 's'] = DictWrapper(
                                bases[0],
                                attr,
                                setter)
                            used_methods.update((attr, setter))
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

        attrs.pop('exclude')
        attrs.pop('preserve')
        return type.__new__(mcs, name, (), attrs)
