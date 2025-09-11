import unittest
from pathlib import Path
import os
import numpy as np

from skycatalogs.skyCatalogs import SkyCatalog, open_catalog
from skycatalogs.utils import Disk
from skycatalogs.objects.base_object import BaseObject


PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent)))
SKYCATALOG_ROOT = os.path.join(PACKAGE_DIR, "skycatalogs", "data", "ci_sample")

SKYCATALOG = os.path.join(SKYCATALOG_ROOT, "trilegal_skyCatalog.yaml")

class TrilegalTester(unittest.TestCase):
    def setUp(self):
        '''
        Open the catalog
        '''
        self._cat = open_catalog(SKYCATALOG, skycatalog_root=SKYCATALOG_ROOT)

    def tearDown(self):
        pass

    def test_disk(self):

        import time


        cat = self._cat

        rgn = Disk(56.6, -36.4, 2200)
        t0 = time.time()
        object_list = cat.get_objects_by_region(rgn, obj_type_set=['trilegal'])
        t_end = time.time()
        print("Time to get disk objects: ", t_end - t0)

        print('Number of collections returned for disk: ',
              object_list.collection_count)
        colls = object_list.get_collections()
        assert(len(object_list.get_collections()) == object_list.collection_count)
        assert(object_list.collection_count == 2)

        # Let it fail if pystellibs isn't available
        # try:
        from pystellibs import BTSettl
        # except ImportError:
        #    return

        from skycatalogs.objects.base_object import load_lsst_bandpasses

        print('Running flux comparison')

        # See if computed fluxes match values in file
        bandpasses = load_lsst_bandpasses()

        for i in range(5):
            obj = object_list[i]
            id = obj.get_native_attribute('id')
            seds = obj.get_observer_sed_components()
            gsobjs = obj.get_gsobject_components()

            # galsim.chromatic.SimpleChromaticTransformation
            chromatic = gsobjs['this_object'] * seds['this_object']

            for b in 'ugrizy':
                calc_flux = chromatic.calculateFlux(bandpasses[b])
                att = f'lsst_flux_{b}'
                cache_flux = obj.get_native_attribute(att)
                assert np.isclose(calc_flux, cache_flux, atol=0) , f'id {id}, band {b}'

if __name__ == '__main__':
    unittest.main()
