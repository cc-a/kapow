from simtk import openmm
# from collections import namedtuple


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

    def __getitem__(self, val):
        try:
            return self.member_wrapper(self.getter(self.parent, val))
        except TypeError:
            return self.member_wrapper(*self.getter(self.parent, val))

    def __setitem__(self, val, args):
        try:
            return self.setter(val, *args)
        except TypeError:
            return self.setter(val, args)

    def __len__(self):
        return self.len_(self.parent)

    def __iter__(self):
        for i in xrange(len(self)):
            yield self[i]
        raise StopIteration

    def _append(self, *val):
        try:
            return self.adder(self.parent, *val)
        except TypeError as e:
            # this can cause uninformative error messages
            print e
            return self.adder(self.parent, val)

    def _pop(self, index):
        val = self[index]
        self.remover(index)
        return val

    def __get__(self, obj, objtype=None):
        self.parent = obj
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

    def __get__(self, obj, objtype=None):
        return self.member_wrapper(property.__get__(self, obj, objtype))

    def __set__(self, obj, value):
        if self.fset is None:
            raise AttributeError("Attribute is read-only")
        try:
            self.fset(obj, *value)
        except TypeError:
            self.fset(obj, value)


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
        return type.__new__(mcs, name, bases, attrs)


class Force(openmm.Force):
    __metaclass__ = Pythonize
    arrays_to_wrap = []
openmm.Force = Force


class NonbondedForce(openmm.NonbondedForce):
    __metaclass__ = Pythonize
    arrays_to_wrap = ['Particle',
                      'Exception']


class System(openmm.System):
    __metaclass__ = Pythonize
    arrays_to_wrap = ['Particle',
                      'Constraint',
                      'Force']
openmm.System = System


class HarmonicBondForce(openmm.HarmonicBondForce, Force):
    __metaclass__ = Pythonize
    arrays_to_wrap = ['Bond']
openmm.HarmonicBondForce = HarmonicBondForce
