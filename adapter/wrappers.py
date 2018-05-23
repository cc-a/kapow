from collections import namedtuple
from functools import update_wrapper
from itertools import chain
import inspect
from .class_map import class_map
from .help_strings import append_help_string, pop_help_string, array_help_string
from .help_strings import dict_help_string


def wrap_method(func):
    """This wrapper provides compatibility by accessing the underlying
    wrapped_object attribute of relevant arguments to provide recognised
    data types in calls to openmm functions.
    """
    def wrapper(*args, **kwargs):
        new_args = []
        for arg in args:
            try:
                new_args.append(arg.wrapped_object)
            except AttributeError:
                new_args.append(arg)
        ret_val = func(*new_args, **kwargs)

        if isinstance(ret_val, tuple):
            # a tuple subtype i.e. a custom member wrapper or tuple
            wrapped_values = [class_map[type(val)].Wrap(val)
                              if type(val) in class_map else val
                              for val in ret_val]
            ret_type = type(ret_val)
            if ret_type is tuple:
                return type(ret_val)(wrapped_values)
            else:
                return type(ret_val)(*wrapped_values)
        else:
            try:
                tmp = class_map[type(ret_val)].Wrap(ret_val)
                return tmp
            except (KeyError, TypeError):
                # ret_val is not a wrapper object
                return ret_val
    return wrapper


class BaseWrapper(object):
    """Fundamental wrapper class.

    Provides default arguments that raise Exceptions if used.
    """
    def default_setter(self, *args, **kwargs):
        raise AttributeError('No setter available for this wrapper')

    def default_adder(self, *args, **kwargs):
        raise AttributeError('No adder available for this wrapper')

    def default_remover(self, *args, **kwargs):
        raise AttributeError('No remover available for this wrapper')

    def default_member_wrapper(self, *args, **kwargs):
        if len(args) == 1:
            ret_val = args[0]
        else:
            ret_val = args

        try:
            if len(ret_val) == 1:
                ret_val = ret_val[0]
        except TypeError:
            pass

        return ret_val


class ArrayWrapper(BaseWrapper):
    """Interface to underlying C++ array mimicing list semantics.

    Provides a pythonic interface to C++ arrays that must be interacted
    with array through getter and setter routines.
    Optional adder and remover routines can be used as well.

    This object provides a list-like interface that wraps the below
    function calls from the corresponding class in openmm:
    Wrapped length function - %s
    Wrapped getters - %s
    Wrapped setters - %s
    Wrapped adder - %s
    Wrapped remover - %s
    """
    def __init__(self, base, len_, getters, setters=[], adder=None,
                 remover=None, member_wrapper=None):
        self.parent = None
        self.len_ = getattr(base, len_)
        self.member_wrapper = member_wrapper

        if adder is None:
            self.adder = self.default_adder
        else:
            self.adder = wrap_method(getattr(base, adder))
            # record the order in which arguments occur in the adder
            # so that we can make sure that the order of getters and
            # setters matches
            spec = inspect.getargspec(getattr(base, adder))
            if spec.defaults is None:
                upper = None
            else:
                upper = -len(spec.defaults)
            args = spec.args[1:upper]  # + spec.kwargs
            # if args == []:
            #     args = ()
            if member_wrapper is None:
                member_name = len_[6:-1]
                if len(args) == 1 and member_name.lower() == args[0].lower():
                    self.member_wrapper = self.default_member_wrapper
                else:
                    self.member_wrapper = namedtuple(member_name, args)

        if len(getters) == 0:
            raise ValueError('At least one getter function must be provided')

        if adder is not None and len(getters) > 1:
            # build getter list such that getters occur in the order
            # matching the adder argument list
            self.getters = []
            for arg in args:
                tmp = [getter for getter in getters
                       if arg.lower() in getter.lower()]
                try:
                    attrname = sorted(tmp, key=lambda x: len(x))[0]
                    self.getters.append(getattr(base, attrname))
                except IndexError:
                    print(tmp, arg, getters)
                    raise
        else:
            self.getters = [getattr(base, getter) for getter in getters]

        if setters == []:
            self.setter = (self.default_setter,)
        else:
            if adder != [] and len(setters) > 1:
                self.setters = [getattr(base, setter) for arg in args
                                for setter in setters if arg in setter.lower()]
            else:
                self.setters = [getattr(base, setter) for setter in setters]

        if remover is None:
            self.remover = self.default_remover
        else:
            self.remover = getattr(base, remover)

        self.member_wrapper = self.default_member_wrapper \
            if self.member_wrapper is None \
            else self.member_wrapper

    @wrap_method
    def __getitem__(self, val):
        if val > len(self) - 1:
            raise IndexError('Array index out of range')

        out = [g(self.parent, val) for g in self.getters]
        for i in range(len(out)):
            if type(out[i]) != list:
                out[i] = [out[i]]
        return self.member_wrapper(*chain(*out))

    def __setitem__(self, val, args):
        try:
            for s in self.setters:
                s(self.parent, val, *args)
            return
        except TypeError:
            pass

        try:
            for s, arg in zip(self.setters, args):
                s(self.parent, val, arg)
            return
        except TypeError:
            pass

        self.setters[0](self.parent, val, args)

    def __len__(self):
        return self.len_(self.parent)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __get__(self, obj, objtype=None):
        try:
            self.parent = obj.wrapped_object
        except AttributeError:
            pass
        return self

    def __set__(self, *args):
        raise AttributeError('Attribute is read-only')

    def __repr__(self):
        try:
            return "<Wrapped C-array containing %d objects>" % len(self)
        except TypeError:
            return "<C-array Wrapper>"


class ValueWrapper(BaseWrapper, property):
    """Wrapper class for a single value or set of values.

    A slight modification of a property. Where values are set with a tuple
    supports argument unpacking for settings expecting multiple values.
    Allows wrapping return values with additional information
    """
    def __init__(self, base, getter, setter='', member_wrapper=None):
        fget = getattr(base, getter)
        if len(inspect.getargspec(fget).args) > 1:
            raise TypeError('Getter function %s must not require an '
                            ' argument other than self' % getter)
        fset = getattr(base, setter, None)
        property.__init__(self, fget=fget, fset=fset,
                          fdel=None, doc=None)
        self.member_wrapper = (self.default_member_wrapper
                               if member_wrapper is None else member_wrapper)

    @wrap_method
    def __get__(self, obj, objtype=None):
        try:
            return self.member_wrapper(
                property.__get__(self, obj, objtype))
        except TypeError:
            return self.member_wrapper(self.fget())
        except AttributeError:
            return self

    def __set__(self, obj, value):
        if self.fset is None:
            raise AttributeError("Attribute is read-only")
        try:
            self.fset(obj.wrapped_object, *value)
        except TypeError:
            self.fset(obj.wrapped_object, value)


class DictWrapper(BaseWrapper, property):
    def __init__(self, fget=None, fset=None, fdel=None, doc=None,
                 frange=None, ffilter=lambda x, y=None: True, fchange=None,
                 fadd=None):
        self.parent = None
        self.member_wrapper = self.default_member_wrapper
        property.__init__(self, fget, fset, fdel, doc)
        self.frange = frange
        self.ffilter = ffilter
        self.fchange = fchange
        self.fadd = fadd

    @wrap_method
    def __getitem__(self, key):
        try:
            return self.member_wrapper(self.fget(self.parent, key))
        except TypeError:
            # static methods
            return self.member_wrapper(self.fget(key))
        except Exception as err:
            raise KeyError(err)

    @wrap_method
    def __setitem__(self, key, val):
        if key not in self.keys() and self.fadd:
            # if key is not present and we have an
            # adder method then use it
            self.fadd(self.parent, key, val)
        else:
            # this may also add new key depending on the setter
            self.fset(self.parent, key, val)

    def setter(self, fset):
        return type(self)(self.fget, fset, self.fdel, self.__doc__,
                          self.frange, self.ffilter, self.fchange, self.fadd)

    def ranger(self, frange):
        return type(self)(self.fget, self.fset, self.fdel, self.__doc__,
                          frange, self.ffilter, self.fchange, self.fadd)

    def keys(self):
        if self.frange and self.ffilter:
            try:
                return [val for val in self.frange(self.parent)
                        if self.ffilter(self.parent, val)]
            except TypeError:
                # staticmethods
                return [val for val in self.frange()
                        if self.ffilter(val)]
        elif self.frange and self.fchange:
            try:
                return [self.fchange(self.parent, val)
                        for val in self.frange(self.parent)]
            except TypeError:
                # staticmethods
                return [self.fchange(val) for val in self.frange()]
        else:
            raise NotImplementedError(
                'cannot do without frange and ffilter yet')

    def iterkeys(self):
        return iter(self.keys())

    def values(self):
        return [self[key] for key in self]

    def itervalues(self):
        return iter(self.values())

    def __iter__(self):
        return self.iterkeys()

    def __len__(self):
        return len(self.keys())

    def __get__(self, obj, objtype=None):
        # Exception block allows staticmethods that don't require
        # a parent object
        try:
            self.parent = obj.wrapped_object
        except AttributeError:
            pass
        return self

    def __set__(self, *args):
        raise AttributeError('Attribute is read-only')


class ClassMethodWrapper(object):
    def __init__(self, cls_method):
        self.cls_method = cls_method

    def __get__(self, obj, objtype=None):
        func = self.cls_method.__get__(obj, objtype)
        wrapped = wrap_method(func)
        update_wrapper(wrapped, func)
        return wrapped


def build_ArrayWrapper(name, base, len_, getters, setters=[], adder=None,
                       remover=None, member_wrapper=None):
    """Create an ArrayWrapper child class with appropriate docstrings
    and methods according to the passed arguments."""
    def append(self, val):
        try:
            return self.adder(self.parent, *val)
        except TypeError:
            # this causes uninformative error messages
            return self.adder(self.parent, val)

    def pop(self, index):
        val = self[index]
        self.remover(self.parent, index)
        return val

    attrs = {}
    if adder is not None:
        adder_method = getattr(base, adder)
        append.__doc__ = append_help_string.format(
            adder_method, adder_method.__doc__)
        attrs['append'] = append

    if remover is not None:
        remover_method = getattr(base, remover)
        pop.__doc__ = pop_help_string.format(
            remover_method, remover_method.__doc__)
        attrs['pop'] = pop

    attrs['__doc__'] = array_help_string.format(
        base, len_, getters, setters, adder, remover)

    cls = type(name, (ArrayWrapper,), attrs)
    return cls(base, len_, getters, setters, adder, remover, member_wrapper)


def build_DictWrapper(name, fget=None, fset=None, fdel=None, doc=None,
                      frange=None, ffilter=lambda x, y=None: True):
    attrs = {'__doc__': dict_help_string.format()}
    cls = type(name, (DictWrapper,), attrs)
    return cls(fget, fset, fdel, doc, frange, ffilter)
