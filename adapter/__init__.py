from simtk import openmm
from simtk.openmm import app
from functools import partial, update_wrapper
# from collections import namedtuple
import re


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
            return args[0]
        else:
            return args


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


class ArrayWrapper(BaseWrapper):
    """Interface to underlying C++ array mimicing list semantics.

    Provides a pythonic interface to C++ arrays that must be interacted
    with array through getter and setter routines.
    Optional adder and remover routines can be used as well.
    """
    def __init__(self, len_, getter, setter=None, adder=None,
                 remover=None, member_wrapper=None):
        self.parent = None
        self.len_ = len_
        self.getter = getter
        self.setter = self.default_setter if setter is None else setter
        self.adder = self.default_adder if adder is None else adder
        self.remover = self.default_remover if remover is None else remover
        self.member_wrapper = (self.default_member_wrapper
                               if member_wrapper is None else member_wrapper)
        if adder:
            self.append = self._append
        if remover:
            self.pop = self._pop

    @wrap_method
    def __getitem__(self, val):
        try:
            return self.member_wrapper(
                self.getter(self.parent, val))
        except TypeError:
            return self.member_wrapper(
                *self.getter(self.parent, val))

    def __setitem__(self, val, args):
        try:
            return self.setter(self.parent, val, *args)
        except TypeError:
            return self.setter(self.parent, val, args)

    def __len__(self):
        return self.len_(self.parent)

    def __iter__(self):
        for i in xrange(len(self)):
            yield self[i]
        raise StopIteration

    @wrap_method
    def _append(self, *val):
        try:
            self.adder(self.parent, *val)
        except TypeError:
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
    def __init__(self, getter, setter=None, member_wrapper=None):
        property.__init__(self, fget=getter, fset=setter,
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


class ForwardFunction(property):
    @wrap_method
    def __get__(self, obj, objtype=None):
        ret_func = partial(self.fget, obj)
        # ret_func.__doc__ = self.fget.__doc__
        update_wrapper(ret_func, self.fget)
        return ret_func

    
class ForwardFunctionAsAttribue(property):
    @wrap_method
    def __get__(self, obj, objtype=None):
        return property.__get__(self, obj, objtype)

    
class Pythonize(type):
    """Metaclass to provide a more Pythonic interface to openmm objects.

    Where values are controlled through getter and setter routines a single
    property interface is created. Where getters, setters and optionally
    add and remove routines are used to interact with array elements a
    more sophisticated interface is created providing List like semantics.
    """
    def __new__(mcs, name, bases, attrs):
        for i in attrs['arrays_to_wrap']:
            method_name = '%ss' % i.lower()
            len_ = getattr(bases[0], 'getNum%ss' % i)
            try:
                getter = next(j for j in bases[0].__dict__
                              if ('get%s' % i) == j)
            except StopIteration:
                try:
                    getter = next(j for j in bases[0].__dict__
                                  if ('get%s' % i) in j)
                except StopIteration:
                    getter = ''
            try:
                setter = next(j for j in bases[0].__dict__
                              if ('set%s' % i) in j)
            except StopIteration:
                setter = ''
            attrs[method_name] = ArrayWrapper(
                len_=len_,
                getter=getattr(bases[0], getter, None),
                setter=getattr(bases[0], setter, None),
                adder=getattr(bases[0], 'add%s' % i))
                # member_wrapper=namedtuple('Particle',
                #                           ('charge','sigma','epsilon')))

        # pull together list of getter methods
        getX = []
        for i in bases[0].__dict__:
            if i.startswith('get'):
                getX.append(i[3:])

        # find getters that have paired setter methods
        pairX = []
        for i in getX:
            try:
                pairX.append(next(j[3:] for j in bases[0].__dict__
                                  if j == 'set%s' % i))
            except StopIteration:
                pass

        # remove elements that are provided as arrays
        for i in pairX[:]:
            for j in attrs['arrays_to_wrap']:
                if j in i:
                    pairX.pop(pairX.index(i))
                    break

        # create attributes for all values
        for i in pairX:
            method_name = i[:1].lower() + i[1:] if i[1:2].islower() else i
            attrs[method_name] = ValueWrapper(getattr(bases[0], 'get%s' % i),
                                              getattr(bases[0], 'set%s' % i))

        attrs.pop('arrays_to_wrap')
        return type.__new__(mcs, name, (), attrs)


# class Force(openmm.Force):
#     __metaclass__ = Pythonize
#     arrays_to_wrap = []
# openmm.Force = Force

# class Force(object):
#     def __init__(self, force):
#         self.force

#     @property
#     def forceGroup(self):
#         return self.force.getForceGroup()


# class NonbondedForce(openmm.NonbondedForce):
#     __metaclass__ = Pythonize
#     arrays_to_wrap = ['Particle',
#                       'Exception']

# class NonbondedForce(object):
#     def __init__(self, nonbondedForce):
#         self.nonbondedForce = nonbondedForce

#     @property
#     def useSwitchingFunction(self):
#         return self.nonbondedForce.getUseSwitchingFunction()

#     @useSwitchingFunction.setter
#     def useSwitchingFunction(self, val):
#         return self.nonbondedForce.setUseSwitchingFunction(val)


# nbf1 = openmm.NonbondedForce()
# nbf2 = NonbondedForce(nbf1)


# def base_init(cls):
#     return init

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

        attrs['__init__'] = init

        attrs_to_add = {}
        for base_class in bases[0].__mro__[:-1]:
            attrs_to_add.update(base_class.__dict__)

        used_methods = []
        get_set_pairs = []
        larger_sets = []
        paired_methods = []
        for method in attrs_to_add:
            if method.startswith('get') and method not in attrs['exclude']:
                words = re.findall('[A-Z][^A-Z]*', method)
                n = ''.join(words[:-1]) if len(words) > 1 else words[0]
                # print "get" + ''.join(words[:-1])
                # has_len = "getNum%ss" % words[0] in attrs_to_add
                # has_setter = ("set%s" % method[3:]) in attrs_to_add
                # has_adder = ("add%s" % words[0]) in attrs_to_add
                # has_remover = ("remove%s" % words[0]) in attrs_to_add
                has_len = "getNum%ss" % n in attrs_to_add
                has_setter = ("set%s" % method[3:]) in attrs_to_add
                has_adder = ("add%s" % n) in attrs_to_add
                has_remover = ("remove%s" % n) in attrs_to_add
                # print method, has_len, has_setter, has_adder
                if has_setter and not (has_len or has_adder or has_remover):
                    if method[4:5].islower():
                        method_name = method[3:4].lower() + method[4:]
                    else:
                        method_name = method[3:]
                    getter = 'get%s' % method[3:]
                    setter = 'set%s' % method[3:]
                    attrs[method_name] = ValueWrapper(
                        getattr(bases[0], getter),
                        getattr(bases[0], setter))
                    used_methods.extend((getter, setter))
                elif has_adder or has_remover:
                    # method_name = words[0][:1].lower() + words[0][1:] if words[0][1:2].islower() else words[0]
                    method_name = n[:1].lower() + n[1:] if n[1:2].islower() else n
                    print method_name, words
                    # len_ = 'getNum%ss' % words[0]
                    # getter = 'get%s' % method[3:]
                    # setter = 'set%s' % method[3:]
                    # adder = 'add%s' % words[0]
                    # remover = 'remove%s' % words[0]
                    len_ = 'getNum%ss' % n
                    getter = 'get%s' % method[3:]
                    setter = 'set%s' % method[3:]
                    adder = 'add%s' % n
                    remover = 'remove%s' % n
                    attrs['%ss' % method_name] = ArrayWrapper(
                        len_=getattr(bases[0], len_),
                        getter=getattr(bases[0], getter),
                        setter=getattr(bases[0], setter, None),
                        adder=getattr(bases[0], adder, None),
                        remover=getattr(bases[0], remover, None))
                    used_methods.extend((len_, getter, setter, adder, remover))
                else:
                    print method
                    

                    
        # if bases[0] == openmm.System:
        #     attrs['virtualSites'] = ArrayWrapper2(getattr(bases[0], 'getVirtualSite'),
        #                                           getattr(bases[0], 'setVirtualSite'))

        # print used_methods
        # for i in bases[0].__dict__:
        for i in attrs_to_add:
            if i not in used_methods and not i.startswith('_') \
               and i not in attrs['exclude']:
                # attrs[i] = ValueWrapper(getattr(bases[0], i))
                attr = getattr(bases[0], i)
                if callable(attr):
                    # if i.startswith('get'):
                    #     if method[4:5].islower():
                    #         method_name = i[3:4].lower() + i[4:]
                    #     else:
                    #         method_name = i[3:]

                    #     attrs[method_name] = ForwardFunctionAsAttribue(attr)
                    # else:
                    attrs[i] = ForwardFunction(attr)
                else:
                    attrs[i] = attr
                used_methods.append(i)


        # print paired_methods, words
        # for method in paired_methods:
        #     words = re.findall('[A-Z][^A-Z]*', method)
        #     print words
        #     sentences = []
        #     for i in xrange(len(words)):
        #         for j in xrange(len(words) - i):
        #             sentences.append(''.join(words[i:i+j+1]))
        #     print method, sentences
        #     tmp = []
        #     for method2 in bases[0].__dict__:
        #         if method not in method2:
        #             tmp.append(longest_common_substring(method, method2))

        #     for lcss in tmp:
        #         if lcss in sentences:
        #             print "matched %s" % lcss

            # for lcss in tmp:
            #     if method.startswith(lcss):
            #         print lcss

        # for i in attrs['arrays_to_wrap']:
        #     len_ = getattr(bases[0], 'getNum%ss' % i)

        #     try:
        #         getter = next(j for j in bases[0].__dict__
        #                       if ('get%s' % i) == j)
        #     except StopIteration:
        #         try:
        #             getter = next(j for j in bases[0].__dict__
        #                           if ('get%s' % i) in j)
        #         except StopIteration:
        #             getter = ''
        #     try:
        #         setter = next(j for j in bases[0].__dict__
        #                       if ('set%s' % i) in j)
        #     except StopIteration:
        #         setter = ''
        #     method_name = '%ss' % i.lower()
        #     attrs[method_name] = ArrayWrapper(
        #         len_=len_,
        #         getter=getattr(bases[0], getter, None),
        #         setter=getattr(bases[0], setter, None),
        #         adder=getattr(bases[0], 'add%s' % i))
        #     used_methods.extend((getter,
        #                          setter,
        #                          'add%s' % i,
        #                          'getNum%ss' % i))

        # # # create attributes for all values
        # # for i in attrs['values_to_wrap']:
        # # pull together list of getter methods
        # getX = []
        # for i in bases[0].__dict__:
        #     if i.startswith('get'):
        #         getX.append(i[3:])

        # # find getters that have paired setter methods
        # pairX = []
        # for i in getX:
        #     try:
        #         pairX.append(next(j[3:] for j in bases[0].__dict__
        #                           if j == 'set%s' % i))
        #     except StopIteration:
        #         pass

        # # remove elements that are provided as arrays
        # for i in pairX[:]:
        #     for j in attrs['arrays_to_wrap']:
        #         if j in i:
        #             pairX.pop(pairX.index(i))
        #             break

        # # create attributes for all values
        # for i in pairX:
        #     method_name = i[:1].lower() + i[1:] if i[1:2].islower() else i
        #     attrs[method_name] = ValueWrapper(getattr(bases[0], 'get%s' % i),
        #                                       getattr(bases[0], 'set%s' % i))

        # for method in bases[0].__dict__:
        #     if method not in used_methods and method not in attrs:
        #         if "__" not in method:
        #             attrs[method] = wrap_method(getattr(bases[0], method))

        # attrs.pop('arrays_to_wrap')
        attrs.pop('exclude')
        return type.__new__(mcs, name, (), attrs)


class NonbondedForce(openmm.NonbondedForce):
    __metaclass__ = Test
    # arrays_to_wrap = ['Particle', 'Exception']
    exclude = []


class System(openmm.System):
    __metaclass__ = Test
    # arrays_to_wrap = ['Particle',
    #                   'Constraint',
    #                   'Force']
    exclude = ['getForces',
               'getVirtualSite',
               'setVirtualSite']


class HarmonicBondForce(openmm.HarmonicBondForce):
    __metaclass__ = Test
    # arrays_to_wrap = ['Bond']
    exclude = []


class CustomNonbondedForce(openmm.CustomNonbondedForce):
    __metaclass__ = Test
    exclude = []

# class HarmonicAngleForce(openmm.HarmonicAngleForce):
#     __metaclass__ = Test
#     # arrays_to_wrap = ['Angle']
#     exclude = []

# class_map = {}
class_map = {openmm.NonbondedForce: NonbondedForce,
             openmm.System: System,
             openmm.HarmonicBondForce: HarmonicBondForce}
# class_map.update({reversed(i) for i in class_map.items()})

prmtop = app.AmberPrmtopFile('../tests/prot_lig1.prmtop')
prmcrd = app.AmberInpcrdFile('../tests/prot_lig1.prmcrd')

residues = list(prmtop.topology.residues())
base = openmm.NonbondedForce()
nbf = NonbondedForce(base)

base = openmm.HarmonicBondForce()
hbf = HarmonicBondForce(base)

base = openmm.System()
sys = System(base)

base = openmm.CustomNonbondedForce('r')
cnbf = CustomNonbondedForce(base)


sys.particles.append(0.)
sys.particles.append(1.)

# vs = openmm.TwoParticleAverageSite(0, 1, 1., 1.)
# sys.virtualSites.append(vs)
