import math
from fractions import Fraction

class squares_sum_solve:
  '''
  here we find integer solutions of q_1 x_1^2 + ... + q_n x_n^2 = c, where q_i and c are rational and positive 
  Input:
  q = [q_1, ... , q_n]
  if offset is given, we look solution for a vector x + offset. Offset is rational.
  '''
  def __init__(self, q, c, offset = None, t=None): 
    self.n = len(q)
    self.offset = [0]*self.n if offset == None else offset
    self.c = c
    self.q = q
    if t == None:
        self.xi_values = self.xi_values_integer
    else:
        self.t = t
        self.xi_values = self.xi_values_quadratic

  def xi_values_integer(self, i, reverse = False): # if not reversed, then partial_sums from 1 to i-1 are used, else from i+1 to n
    start = -math.sqrt(self.c/self.q[i])-self.offset[i]
    stop = math.sqrt(self.c/self.q[i])-self.offset[i]
    start, stop = int(math.floor(start)), int(math.floor(stop))
    return {xi: self.q[i]*(xi+self.offset[i])**2 for xi in range(start, stop+2)}

  def xi_values_quadratic(self, i):
            Cq = self.c/self.q[i]
            #print 'Cq =', float(Aut(Cq))
            start_a = math.floor((-math.sqrt(Cq) - math.sqrt(abs(Aut(Cq))))/2)
            stop_a = math.floor((math.sqrt(Cq) + math.sqrt(abs(Aut(Cq))))/2)
            start_a, stop_a = int(math.floor(start_a)), int(math.floor(stop_a))
            Ai = range(start_a - 1, stop_a + 2 )
            xi_pairs = ((ai, bi) for ai in Ai for bi in range(int(math.floor((-math.sqrt(Cq) - ai)/t)), int(math.floor((math.sqrt(Cq) - ai)/t))+2))
            return {x : self.q[i]*(x[0] + x[1]*self.t)**2 for x in xi_pairs}  # possible values of qi*xi^2

  #print('squares_sum_solve, sums:',sums[n-1])

  def is_possible_partial_sum(self, s):
    return s <= self.c


  def compute_partial_sums(self):
    self.partial_sums = [self.xi_values(0).values()]
    for i in range(1,self.n):
      summands = self.xi_values(i).values()
      sums = (p+q for p in self.partial_sums[-1] for q in summands if self.is_possible_partial_sum(p+q)) # sums of all pairs
      sums = list(dict.fromkeys(sums)) # removed duplicates
      self.partial_sums.append(sums) 

  def run(self):
    self.compute_partial_sums()
    if self.c not in self.partial_sums[self.n-1]:
      print('squares_sum_solve: no solutions found')
      return
    for s in self.solutions(self.n-1,[0]*self.n,self.c):
        yield s
    
    

  def solutions(self, i, x, remainder):
    if i==0:
      x_1_values = [k for k, v in self.xi_values(0).items() if v == remainder]
      if len(x_1_values)>2: # sanity check
        print('ERROR (internal): wrong sum\n','q: ',self.q,'\nc: ',self.c,'\noffset: ',self.offset, '\nremainder: ',remainder, '\nx1: ',x_1_values)
      for x_1 in x_1_values:
        x[0] = x_1
        yield x
    else:
      for xi, si in self.xi_values(i).items(): # si = qi*xi*xi
        if remainder - si in self.partial_sums[i-1]:
          x[i]=xi
          for s in self.solutions(i-1, x.copy(), remainder - si):
                yield s

  
if __name__ == '__main__':
  q = [1,Fraction(1,2),3,4]
  offset = [0,0,3,4]
  print([s for s in squares_sum_solve(q, 2).run()])
  print([s for s in squares_sum_solve(q, Fraction(1,2), offset = offset).run()])