from simtk import openmm
# from simtk.openmm import app
from functools import partial
from six import with_metaclass
from collections import namedtuple, defaultdict
import inspect
import sys
from .pythonize import Pythonize
from .wrappers import DictWrapper, ValueWrapper, build_ArrayWrapper
from .class_map import class_map
from . import app

# TODO

# test reload logic

# possible for getters with a getNum and a getXXXName could build a dictwrapper

# what is going on with the integrators?


def print_args(func):
    """For debugging"""
    def wrapped(*args):
        print(args)
        return func(*args)
    return wrapped


class Context(with_metaclass(Pythonize, openmm.Context)):
    exclude = ['getParameter', 'setParameter', 'getParameters']
    preserve = []

    @DictWrapper
    def parameters(self, key):
        return openmm.Context.getParameter(self, key)

    @parameters.setter
    def parameters(self, key, val):
        return openmm.Context.setParameter(self, key, val)

    @parameters.ranger
    def parameters(self):
        return openmm.Context.getParameters(self)


class System(with_metaclass(Pythonize, openmm.System)):
    exclude = ['getForces', 'getVirtualSite',
               'setVirtualSite', 'isVirtualSite']
    preserve = []

    @partial(DictWrapper, ffilter=openmm.System.isVirtualSite,
             frange=lambda x: range(openmm.System.getNumParticles(x)))
    def virtualSites(self, key):
        return openmm.System.getVirtualSite(self, key)

    @virtualSites.setter
    def virtualSites(self, key, value):
        return openmm.System.setVirtualSite(self, key, value)


class Platform(with_metaclass(Pythonize, openmm.Platform)):
    exclude = ['getPropertyDefaultValue', 'setPropertyDefaultValue',
               'getNumPlatforms', 'getPlatform', 'getPlatformByName']
    preserve = ['getPropertyValue', 'setPropertyValue']

    @DictWrapper
    def propertyDefaults(self, key):
        return openmm.Platform.getPropertyDefaultValue(self, key)

    @propertyDefaults.setter
    def propertyDefaults(self, key, value):
        return openmm.Platform.setPropertyDefaultValue(self, key, value)

    propertyDefaults = propertyDefaults.ranger(
        openmm.Platform.getPropertyNames)

    @DictWrapper
    def platformsByName(key):
        return openmm.Platform.getPlatformByName(key)

    @platformsByName.ranger
    def platformsByName():
        return [openmm.Platform.getPlatform(i).getName()
                for i in range(openmm.Platform.getNumPlatforms())]


class AmoebaMultipoleForce(with_metaclass(Pythonize, openmm.AmoebaMultipoleForce)):
    exclude = ['getCovalentMap', 'setCovalentMap', 'getCovalentMaps']
    preserve = []

    @DictWrapper
    def covalentMaps(self, key):
        return openmm.AmoebaMultipoleForce.getCovalentMaps(self, key)

    @covalentMaps.setter
    def covalentMaps(self, key, value):
        if len(value) != 8:
            raise TypeError('Covalent Map must be a tuple of length 8')
        for i, val in enumerate(value):
            openmm.AmoebaMultipoleForce.setCovalentMap(self, key, i, val)

    @covalentMaps.ranger
    def covalentMaps(self):
        return range(openmm.AmoebaMultipoleForce.getNumMultipoles(self))


class TwoParticleAverageSite(with_metaclass(Pythonize, openmm.TwoParticleAverageSite)):
    exclude = ['getNumParticles', 'getParticle', 'getWeight']
    preserve = []
    particles = build_ArrayWrapper(
        'particles',
        openmm.TwoParticleAverageSite,
        'getNumParticles',
        ['getParticle', 'getWeight'],
        member_wrapper=namedtuple('Particle', ('index', 'weight')))


class CustomCentroidBondForce(with_metaclass(Pythonize, openmm.CustomCentroidBondForce)):
    exclude = ['getNumGroupsPerBond']
    preserve = []
    numGroupsPerBond = ValueWrapper(
        openmm.CustomCentroidBondForce, 'getNumGroupsPerBond')


class CustomCompoundBondForce(with_metaclass(Pythonize, openmm.CustomCompoundBondForce)):
    exclude = ['getNumParticlesPerBond']
    preserve = []
    numParticlesPerBond = ValueWrapper(
        openmm.CustomCompoundBondForce, 'getNumParticlesPerBond')


class CustomManyParticleForce(with_metaclass(Pythonize, openmm.CustomManyParticleForce)):
    exclude = ['getNumParticlesPerSet']
    preserve = []
    numParticlesPerSet = ValueWrapper(
        openmm.CustomManyParticleForce, 'getNumParticlesPerSet')


exclusions = defaultdict(
    lambda: [],
    {'System': ['getForces'],  # convenience function that is not needed
     'CustomManyParticleForce': [
         'getTypeFilter'],  # index, value pair
     'Platform': ['getPropertyValue'] # should this be a preserve for get and set?
    } 
)

skip = [
    'RPMDIntegrator',
    'CustomIntegrator',
    'DualAMDIntegrator',
    'MTSIntegrator',
    'AMDForceGroupIntegrator',
    'AMDIntegrator',
    'AmoebaTorsionTorsionForce',
    'AmoebaVdwForce'
]

module = sys.modules[__name__]
for name in [n for n in dir(openmm) if inspect.isclass(getattr(openmm, n))]:
    if name in skip:
        continue
    # does this logic cause problems for reloads?
    if getattr(module, name, None) is None:
        openmm_class = getattr(openmm, name)
        new_class = Pythonize(name,
                              [openmm_class],
                              {'exclude': exclusions[name], 'preserve': []})
        setattr(module, name, new_class)
    class_map[getattr(openmm, name)] = getattr(module, name)
