from vinal import *
import unittest
import yaml


def block(lattice_type, n=0):
    if lattice_type == 'U':
        return [[0,1],[1,0]]
    if lattice_type == 'A':
        return [[max(2-abs(i-j),0)  for j in range(n)] for i in range(n)]

class VinalTestCase(unittest.TestCase):
    def setUp(self):
      with open('lattices.yaml', "r") as f:
        self.lattices = yaml.load(f)['lattices']

    def test_En(self):
      for L in self.lattices:
        M = blocks(L['blocks'])
        self.assertEqual(VinAl(M).En, L['En'],
                        f'incorrect En for entry {L}')
        self.assertTrue(VinAl(L['Q']).is_root(L['root_example']), f'incorrect root_example for entry {L}')


if __name__ == '__main__':
    unittest.main()


# по идее, можно запускать в терминале на компе командой python -m unittest tests, но у меня не получилось

# стоит добавить все тесты из check_validity в sagemath-версии, а также из init