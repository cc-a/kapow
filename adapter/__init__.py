from simtk import openmm
from simtk.openmm import app
from functools import partial, update_wrapper
from collections import namedtuple, defaultdict
import inspect
import sys

# TODO

# Virtual site getWeight function

# DictWrapper dynamic key determination


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

        # try:
        #     args = (args[0].wrapped_object,) + args[1:]
        # except AttributeError:
        #     pass
        # ret_val = func(*args, **kwargs)
        try:
            return class_map[ret_val.__class__](ret_val)
        except (KeyError, TypeError):
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


# class ArrayWrapper2(BaseWrapper):
#     def __init__(self, getter, adder):
#         self.parent = None
#         self.values = []
#         self.getter = getter
#         self.adder = adder

#     def __get__(self, obj, objtype=None):
#         try:
#             self.parent = obj.wrapped_object
#         except AttributeError:
#             pass
#         return self

#     def __getitem__(self, val):
#         return self.values[val]

#     def append(self, *val):
#         self.values.append(val)
#         self._append(self, *val)

#     @wrap_method
#     def _append(self, *val):
#         try:
#             self.adder(self.parent, *val)
#         except TypeError:
#             # this causes uninformative error messages
#             self.adder(self.parent, val)


# class ArrayWrapper(BaseWrapper):
#     """Interface to underlying C++ array mimicing list semantics.

#     Provides a pythonic interface to C++ arrays that must be interacted
#     with array through getter and setter routines.
#     Optional adder and remover routines can be used as well.
#     """
#     def __init__(self, base, len_, getters, setters=None, adder=None,
#                  remover=None, member_wrapper=None):
#         self.parent = None
#         self.len_ = getattr(base, len_)
#         try:
#             len(getter)
#             self.getter = getter
#         except TypeError:
#             self.getter = (getter,)

#         if setter is None:
#             self.setter = self.default_setter
#         else:
#             try:
#                 len(setter)
#                 self.setter = setter
#             except TypeError:
#                 self.setter = (setter,)
#         # self.setter = self.default_setter if setter is None else setter
#         self.adder = self.default_adder if adder is None else adder
#         self.remover = self.default_remover if remover is None else remover
#         self.member_wrapper = (self.default_member_wrapper
#                                if member_wrapper is None else member_wrapper)
#         if adder:
#             self.append = self._append
#         if remover:
#             self.pop = self._pop

#     @wrap_method
#     def __getitem__(self, val):
#         try:
#             # return self.member_wrapper(
#             #     self.getter(self.parent, val))
#             return self.member_wrapper([g(self.parent, val)
#                                         for g in self.getter])
#         except TypeError:
#             return self.member_wrapper(
#                 *self.getter(self.parent, val))

#     def __setitem__(self, val, args):
#         try:
#             # return self.setter(self.parent, val, *args)
#             for s in self.setter:
#                 s(self.parent, val, *args)
#             return
#         except TypeError:
#             pass

#         try:
#             for s, arg in zip(self.setter, args):
#                 s(self.parent, val, arg)
#             return
#         except TypeError:
#             pass

#         self.setter[0](self.parent, val, args)

#     def __len__(self):
#         return self.len_(self.parent)

#     def __iter__(self):
#         for i in xrange(len(self)):
#             yield self[i]
#         raise StopIteration

#     @wrap_method
#     def _append(self, *val):
#         try:
#             self.adder(self.parent, *val)
#         except TypeError:
#             # this causes uninformative error messages
#             self.adder(self.parent, val)

#     def _pop(self, index):
#         val = self[index]
#         self.remover(self.parent, index)
#         return val

#     def __get__(self, obj, objtype=None):
#         try:
#             self.parent = obj.wrapped_object
#         except AttributeError:
#             pass
#         return self

#     def __set__(self, *args):
#         raise AttributeError('Attribute is read-only')

#     def __repr__(self):
#         try:
#             return "<Wrapped C-array containing %d objects>" % len(self)
#         except TypeError:
#             return "<C-array Wrapper>"

def print_args(func):
    def wrapped(*args):
        print args
        return func(*args)
    return wrapped


class ArrayWrapper(BaseWrapper):
    """Interface to underlying C++ array mimicing list semantics.

    Provides a pythonic interface to C++ arrays that must be interacted
    with array through getter and setter routines.
    Optional adder and remover routines can be used as well.
    """
    def __init__(self, base, len_, getters, setters=None, adder=None,
                 remover=None, member_wrapper=None):
        self.parent = None
        self.len_ = getattr(base, len_)
        self.member_wrapper = None

        if adder == []:
            self.adder = self.default_adder
        else:
            self.adder = getattr(base, adder[0])
            # record the order in which arguments occur in the adder
            # so that we can make sure that the order of getters and
            # setters matches
            spec= inspect.getargspec(self.adder)
            args = spec.args[1:] #+ spec.kwargs
            # if args == []:
            #     args = ()
            # self.member_wrapper = namedtuple(len_[6:-1], args)

        if adder != [] and len(getters) > 1:
            self.getters = [getattr(base, getter) for arg in args
                            for getter in getters
                            if arg.lower() in getter.lower()]
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

        if remover == []:
            self.remover = self.default_remover
        else:
            self.remover = getattr(base, remover[0])

        self.member_wrapper = (self.default_member_wrapper
                               if self.member_wrapper is None else self.member_wrapper)
        # self.member_wrapper = self.member_wrapper or self.default_member_wrapper
        if adder is not None:
            self.append = self._append
            # update_wrapper(self.append, self.adder, assigned=('__doc__',))
        if remover is not None:
            self.pop = self._pop

    @wrap_method
    def __getitem__(self, val):
        # try:
            # return self.member_wrapper(
            #     self.getter(self.parent, val))
        # print [g(self.parent, val)
        #        for g in self.getters]
        
        return self.member_wrapper(tuple(g(self.parent, val)
                                         for g in self.getters))
        # except TypeError:
        #     return self.member_wrapper(
        #         *self.getter(self.parent, val))

    def __setitem__(self, val, args):
        try:
            # return self.setter(self.parent, val, *args)
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
        for i in xrange(len(self)):
            yield self[i]
        raise StopIteration

    @wrap_method
    def _append(self, val):
        try:
            self.adder(self.parent, *val)
        except TypeError as e:
            # print (e)
            # this causes uninformative error messages
            self.adder(self.parent, val)

    def _pop(self, index):
        val = self[index]
        self.remover(self.parent, index)
        return val

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
    """At the current time this is fundamentally flawed. Adding members
    directly to the wrapped object will not update self._keys.
    Ultimately the key list must be built dynamically from the wrapped object.
    """
    def __init__(self, *args, **kwargs):
        self.parent = None
        self.member_wrapper = self.default_member_wrapper
        property.__init__(self, *args)
        self._keys = []

    #     self.getter = getattr(base, getter)
    #     self.setter = getattr(base, setter)
    # @staticmethod
    # def withFilter(filter_func):
    #     return partial(DictWrapper, filter_func=filter_func)

    @wrap_method
    def __getitem__(self, key):
        # try:
            # return self.member_wrapper(
            #     self.getter(self.parent, val))
        # print [g(self.parent, val)
        #        for g in self.getters]
        
        # return self.member_wrapper(self.getter(self.parent, val))
        try:
            return self.member_wrapper(self.fget(self.parent, key))
        except Exception as err:
            raise KeyError(err)
        # except TypeError:
        #     return self.member_wrapper(
        #         *self.getter(self.parent, val))

    @wrap_method
    def __setitem__(self, key, val):
        self.fset(self.parent, key, val)
        self._keys.append(key)

    def keys(self):
        return list(self._keys)

    def iterkeys(self):
        return iter(self._keys)

    def values(self):
        return [self[key] for key in self]

    def itervalues(self):
        return iter(self.values())

    def __iter__(self):
        return self.iterkeys()

    def __len__(self):
        return len(self.keys())

    def __get__(self, obj, objtype=None):
        try:
            self.parent = obj.wrapped_object
        except AttributeError:
            pass
        return self

    def __set__(self, *args):
        raise AttributeError('Attribute is read-only')

def init(self, inst):
    self.wrapped_object = inst


# taken from wikipedia - so yeah...
def longest_common_substring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


class Test(type):
    """Metaclass to provide a more Pythonic interface to openmm objects.

    Where values are controlled through getter and setter routines a single
    property interface is created. Where getters, setters and optionally
    add and remove routines are used to interact with array elements a
    more sophisticated interface is created providing List like semantics.
    """
    def __new__(mcs, name, bases, attrs):

        if '__init__' not in attrs:
            attrs['__init__'] = init
        attrs_to_add = {}
        for base_class in bases[0].__mro__[:-1]:
            attrs_to_add.update(base_class.__dict__)

        for method in attrs['exclude']:
            attrs_to_add.pop(method)

        all_array_roots = [attr[6:-1]
                           for attr in attrs_to_add.keys()
                           if 'getNum' in attr]

        array_getters = {}
        for root in all_array_roots:
            array_getters[root] = [attr for attr in attrs_to_add.keys()
                                   if 'get' + root in attr]

        array_setters = {}
        for root in all_array_roots:
            array_setters[root] = [attr for attr in attrs_to_add.keys()
                                   if 'set' + root in attr]

        array_adders = {}
        for root in all_array_roots:
            # try:
            array_adders[root] = [attr for attr in attrs_to_add.keys()
                                  if 'add' + root in attr]
            # except IndexError:
            #     array_adders[root] = None

        array_removers = {}
        for root in all_array_roots:
            # try:
            array_removers[root] = [attr for attr in attrs_to_add.keys()
                                    if 'remove' + root in attr]
            # except IndexError:
            #     array_removers[root] = None

        for root in all_array_roots:
            method_name = root[0].lower() + root[1:] + 's'
            attrs[method_name] = ArrayWrapper(
                bases[0],
                len_='getNum%ss' % root,
                getters=array_getters[root],
                setters=array_setters[root],
                adder=array_adders[root],
                remover=array_removers[root])
        used_methods = set(method for method_type in (array_adders,
                                                      array_removers,
                                                      array_getters,
                                                      array_setters)
                           for methods in method_type.values()
                           for method in methods)
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
        return type.__new__(mcs, name, (), attrs)


class System(openmm.System):
    __metaclass__ = Test
    exclude = ['getForces', 'getVirtualSite', 'setVirtualSite', 'isVirtualSite']
    def __init__(self, inst):
        self.wrapped_object = inst

    @DictWrapper
    def virtualSites(self, key):
        return openmm.System.getVirtualSite(self, key)

    @virtualSites.setter
    def virtualSites(self, key, value):
        return openmm.System.setVirtualSite(self, key, value)

exclusions = defaultdict(
    lambda: [],
    {'System': ['getForces'],  # convenience function that is not needed
                # 'getVirtualSite',
                # 'setVirtualSite'],  # index, value pair
     'AmoebaMultipoleForce': ['getCovalentMap', # (index, typeid), value pair
                              'setCovalentMap'], # also getCovalentMaps
     'Context': ['getParameter', 'setParameter', # name, value pair, also getParameters
                 'getState'],
     'CustomCentroidBondForce': ['getNumGroupsPerBond'], # name mangling
     'CustomCompoundBondForce': ['getNumParticlesPerBond'], # name mangling
     'CustomManyParticleForce': ['getNumParticlesPerSet', # name mangling
                                 'getTypeFilter'], # index, value pair
     'Platform': ['getPropertyDefaultValue', 'getPropertyValue']}) # name, value pair

module = sys.modules[__name__]
class_map = {}
for name in [n for n in dir(openmm) if inspect.isclass(getattr(openmm, n))]:
    # does this logic cause problems for reloads?
    if getattr(module, name, None) is None:
        openmm_class = getattr(openmm, name)
        new_class = Test(name,
                         [openmm_class],
                         {'exclude': exclusions[name]})
        setattr(module, name, new_class)
        class_map[openmm_class] = new_class
    else:
        print 'skipping', name


# class CustomNonbondedForce(openmm.CustomNonbondedForce):
#     __metaclass__ = Test
#     exclude = ()
#     # tabulatedFunctions = ArrayWrapper(
#     #     len_=openmm.CustomNonbondedForce.getNumFunctions
#     #     getter=(openmm.CustomNonbondedForce.getTabulatedFunction,
#     #             openmm.CustomNonbondedForce.getFunctionParameters,
#     #             openmm.CustomNonbondedForce.getTabulatedFunctionName),
#     #     setter=(openmm.CustomNonbondedForce.setF))


prmtop = app.AmberPrmtopFile('../tests/prot_lig1.prmtop')
prmcrd = app.AmberInpcrdFile('../tests/prot_lig1.prmcrd')

from simtk.unit import picoseconds
system = prmtop.createSystem()
vs = openmm.TwoParticleAverageSite(0, 1, 0.5, 0.5)
system.setVirtualSite(0, vs)
print id(vs), id(system.getVirtualSite(0))

# wrapped = System(system)
# integrator = openmm.VerletIntegrator(0.001 * picoseconds)
# context = openmm.Context(system, integrator)


# residues = list(prmtop.topology.residues())

# base = openmm.NonbondedForce()
# nbf = NonbondedForce(base)

# base = openmm.CustomNonbondedForce('r')
# cnbf = CustomNonbondedForce(base)

# tfunc = openmm.Continuous1DFunction((0,1), 0., 1.)

# base.addTabulatedFunction('first', tfunc)
# base.addFunction('secend', (0, 1), 0, 1)

# cnbf.globalParameters.append(('jon', 0))
# print cnbf.globalParameters[0]
# print len(cnbf.globalParameters), cnbf.globalParameters[0]
# cnbf.globalParameters[0] = ('chris', 1)
# print len(cnbf.globalParameters), cnbf.globalParameters[0]

# cnbf.perParticleParameters.append(('jon',))
# cnbf.particles.append(())
# print cnbf.particles[0]
# cnbf.particles.append(((1,),))
# print cnbf.particles[1]
# # hbf = HarmonicBondForce(base)

# # base = openmm.System()
# # sys = System(base)

# base = openmm.CustomNonbondedForce('r')
# cnbf = CustomNonbondedForce(base)

# tfunc = openmm.Continuous1DFunction((0,1), 0., 1.)

# base.addTabulatedFunction('first', tfunc)
# base.addFunction('secend', (0, 1), 0, 1)

# cnbf.globalParameters.append(('jon', 0))
# print cnbf.globalParameters[0]
# print len(cnbf.globalParameters), cnbf.globalParameters[0]
# cnbf.globalParameters[0] = ('chris', 1)
# print len(cnbf.globalParameters), cnbf.globalParameters[0]

# cnbf.perParticleParameters.append(('jon',))
# cnbf.particles.append(())
# print cnbf.particles[0]
# cnbf.particles.append(([1],))
# print cnbf.particles[1]


