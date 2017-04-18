import unittest
import adapter
from simtk import openmm
from simtk.unit import nanometer, elementary_charge, kilojoule_per_mole

class TestWrappedClasses(unittest.TestCase):
    def testNonbondedForce(self):
        inst = adapter.NonbondedForce()
        inst.reactionFieldDielectric = 1.0
        self.assertEqual(inst.reactionFieldDielectric, 1.0)
        
        inst.useSwitchingFunction
        inst.useSwitchingFunction = False
        self.assertEqual(inst.useSwitchingFunction, False)
        inst.useSwitchingFunction = True
        self.assertEqual(inst.useSwitchingFunction, True)

        inst.PMEParameters = 1.0 / nanometer, 1, 1, 1
        self.assertEqual(inst.PMEParameters,
                         [1.0 / nanometer, 1, 1, 1])

        inst.particles.append(0.0, 0.0, 0.0)
        self.assertEqual(inst.particles[0],
                         [0. * elementary_charge,
                          0. * nanometer,
                          0. * kilojoule_per_mole])
        self.assertEqual(inst.particles[0],
                         inst.getParticleParameters(0))
        inst.exceptions.append(1,2, 0.0, 0.0, 0.0)
        self.assertEqual(inst.exceptions[0],
                         [1, 2,
                          0. * elementary_charge**2,
                          0. * nanometer,
                          0. * kilojoule_per_mole])
        self.assertEqual(inst.exceptions[0], inst.getExceptionParameters(0))

        # forceGroup is inherited from the Force class
        inst.forceGroup = 2
        self.assertEqual(inst.forceGroup, 2)

    def testHarmonicBondForce(self):
        inst = adapter.HarmonicBondForce()
        inst.bonds.append(0, 1, 0. * nanometer,
                          0. * kilojoule_per_mole / nanometer**2)
        self.assertEqual(inst.bonds[0],
                         [0, 1, 0. * nanometer,
                          0. * kilojoule_per_mole / nanometer**2])
        self.assertEqual(inst.bonds[0], inst.getBondParameters(0))

        # forceGroup is inherited from the Force class
        inst.forceGroup = 2
        self.assertEqual(inst.forceGroup, 2)


# prmtop = openmm.app.AmberPrmtopFile('../multisim/plugin/tests/prot_lig1.prmtop')
# prmcrd = openmm.app.AmberInpcrdFile('../multisim/plugin/tests/prot_lig1.prmcrd')

# residues = list(prmtop.topology.residues())
# system = prmtop.createSystem(constraints=openmm.app.HBonds,
#                              removeCMMotion=True)

# #particles # need virtual site tests
# wrapped.particles
# wrapped.particles[0]
# wrapped.particles[1] = 0.0
# wrapped.particles.append(1.0)
# vs = openmm.TwoParticleAverageSite(0, 1, 0.5, 0.5)
# wrapped.particles.append((1.0, vs))
# for p in wrapped.particles: pass

if __name__ == '__main__':
    unittest.main()



    
