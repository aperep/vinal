import vinal
import unittest
import yaml

class VinalTestCase(unittest.TestCase):
    def setUp(self):
      with open('lattices.yaml', "r") as f:
        self.lattices = yaml.load(f)['lattices']

    def test_En(self):
      for L in self.lattices:
        self.assertEqual(vinal.VinAl(L['Q']).En, L['En'],
                         'incorrect En')

if __name__ == '__main__':
    unittest.main()


# по идее, можно запускать в терминале на компе командой python -m unittest tests, но у меня не получилось

# стоит добавить все тесты из check_validity в sagemath-версии, а также из init