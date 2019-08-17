import math
from fractions import Fraction
# 
def squares_sum_solve(q, c, offset = None): 
  '''
  here we find integer solutions of q_1 x_1^2 + ... + q_n x_n^2 = c, where q_i and c are rational and positive 
  Input:
  q = [q_1, ... , q_n]
  if offset is given, we look solution for a vector x + offset. Offset is rational.
  '''
  # 
  n = len(q)
  if offset == None:
    offset = [0]*n

  def possible_x_range(i):
    return range( math.floor(-math.sqrt(c/q[i])-offset[i]), math.floor(math.sqrt(c/q[i])-offset[i])+1 )
  summands = [{x:q[i]*(x+offset[i])**2 for x in possible_x_range(i)} for i in range(n)] # lists of possible values of qi*xi^2
  sums = [set(summands[0].values())] # ith element consists of all possible sums of first i summands
  for i in range(1,n):
    sums.append({p+q for p in sums[i-1] for q in summands[i].values() if p+q <=c})
  #print('squares_sum_solve, sums:',sums[n-1])
  if c not in sums[n-1]:
    #print('squares_sum_solve: no solutions found')
    pass

  def solutions(i, x, remainder):
    if i==0:
      x_1_values = [k for k, v in summands[0].items() if v == remainder]
      if len(x_1_values)>2: # sanity check
        print('ERROR (internal): wrong sum\n','q: ',q,'\nc: ',c,'\noffset: ',offset,
          '\nsummands[0]: ',summands[0], '\nremainder: ',remainder, '\nx1: ',x_1_values)
      for x_1 in x_1_values:
        x[0] = x_1
        yield x
    else:
      for xi, si in summands[i].items(): # si = qi*xi*xi
        if remainder - si in sums[i-1]:
          x[i]=xi
          yield from solutions(i-1, x.copy(), remainder - si)

  yield from solutions(n-1,[0]*n,c)
  

if __name__ == '__main__':
  q = [1,Fraction(1,2),3,4]
  offset = [0,0,3,4]
  print([s for s in squares_sum_solve(q, 2)])
  print([s for s in squares_sum_solve(q, Fraction(1,2), offset = offset)])