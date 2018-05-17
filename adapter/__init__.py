from simtk import openmm
# from simtk.openmm import app
from functools import partial, update_wrapper
from collections import namedtuple, defaultdict
import inspect
import sys
from .pythonize import Pythonize
from .wrappers import ArrayWrapper, DictWrapper, ValueWrapper
from class_map import class_map

# TODO

# test reload logic

# possible for getters with a getNum and a getXXXName could build a dictwrapper


def print_args(func):
    """For debugging"""
    def wrapped(*args):
        print args
        return func(*args)
    return wrapped

# # taken from wikipedia - so yeah...
# def longest_common_substring(s1, s2):
#     m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
#     longest, x_longest = 0, 0
#     for x in xrange(1, 1 + len(s1)):
#         for y in xrange(1, 1 + len(s2)):
#             if s1[x - 1] == s2[y - 1]:
#                 m[x][y] = m[x - 1][y - 1] + 1
#                 if m[x][y] > longest:
#                     longest = m[x][y]
#                     x_longest = x
#             else:
#                 m[x][y] = 0
#     return s1[x_longest - longest: x_longest]


class System(openmm.System):
    __metaclass__ = Pythonize
    exclude = ['getForces', 'getVirtualSite',
               'setVirtualSite', 'isVirtualSite']
    preserve = []

    @partial(DictWrapper, ffilter=openmm.System.isVirtualSite,
             frange=lambda x: xrange(openmm.System.getNumParticles(x)))
    def virtualSites(self, key):
        return openmm.System.getVirtualSite(self, key)

    @virtualSites.setter
    def virtualSites(self, key, value):
        return openmm.System.setVirtualSite(self, key, value)


class Platform(openmm.Platform):
    __metaclass__ = Pythonize
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
                for i in xrange(openmm.Platform.getNumPlatforms())]


class AmoebaMultipoleForce(openmm.AmoebaMultipoleForce):
    __metaclass__ = Pythonize
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
        return xrange(openmm.AmoebaMultipoleForce.getNumMultipoles(self))


class TwoParticleAverageSite(openmm.TwoParticleAverageSite):
    __metaclass__ = Pythonize
    exclude = ['getNumParticles', 'getParticle', 'getWeight']
    preserve = []
    particles = ArrayWrapper(
        openmm.TwoParticleAverageSite,
        'getNumParticles',
        ['getParticle', 'getWeight'])


exclusions = defaultdict(
    lambda: [],
    {'System': ['getForces'],  # convenience function that is not needed
     'CustomCentroidBondForce': ['getNumGroupsPerBond'],  # name mangling
     'CustomCompoundBondForce': ['getNumParticlesPerBond'],  # name mangling
     'CustomManyParticleForce': ['getNumParticlesPerSet',  # name mangling
                                 'getTypeFilter'],  # index, value pair
     'Platform': ['getPropertyValue']})

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
    else:
        print 'skipping', name
    class_map[getattr(openmm, name)] = getattr(module, name)


# class CustomNonbondedForce(openmm.CustomNonbondedForce):
#     __metaclass__ = Pythonize
#     exclude = ()
#     # tabulatedFunctions = ArrayWrapper(
#     #     len_=openmm.CustomNonbondedForce.getNumFunctions
#     #     getter=(openmm.CustomNonbondedForce.getTabulatedFunction,
#     #             openmm.CustomNonbondedForce.getFunctionParameters,
#     #             openmm.CustomNonbondedForce.getTabulatedFunctionName),
#     #     setter=(openmm.CustomNonbondedForce.setF))


# prmtop = app.AmberPrmtopFile('../tests/prot_lig1.prmtop')
# prmcrd = app.AmberInpcrdFile('../tests/prot_lig1.prmcrd')

# from simtk.unit import picoseconds
# system = prmtop.createSystem()
# integrator = openmm.VerletIntegrator(0.001 * picoseconds)
# context = openmm.Context(system, integrator)
