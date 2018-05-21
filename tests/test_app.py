from adapter import app as app
from simtk.openmm import app as oapp
import unittest


class Test(unittest.TestCase):
    def testAttributes(self):
        """A very simple test that makes sure there is a 1-1 mapping
        between the attributes of the original and wrapped modules.
        """
        for name in dir(oapp):
            if name in app.skip or name.startswith('_'):
                continue
            obj = getattr(app, name)
            oobj = getattr(oapp, name)

            for attr in dir(oobj):
                if not attr.startswith('_'):
                    try:
                        getattr(obj, attr)
                    except AttributeError:
                        print name, obj, attr
                        raise
