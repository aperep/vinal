import math

# 
def squares_sum_solve(q, c, offset = None): 
  '''
  here we find integer solutions of q_1^2 x_1 + ... + q_n^2 x_n = c, where q_i and c are integer and positive
  Input:
  q = [q_1, ... , q_n]
  if offset is given, we look solution for a vector q + offset. Offset can be rational. Not implemented.
  '''
  # 
  n = len(q)
  if offset == None:
    offset = [0]*n

  def possible_x_range(i):
    return range(math.floor(-math.sqrt(c/q[i])-offset[i]),math.floor(math.sqrt(c/q[i])-offset[i])+1)
  summands = [{x:q[i]*(x+offset[i])**2 for x in possible_x_range(i)} for i in range(n)] # lists of possible values of ai*xi^2
  sums = [set(summands[0].values())]
  for i in range(1,n):
    sums.append({p+q for p in sums[i-1] for q in summands[i].values() if p+q <=c})
  print('sums:',sums[n-1])
  if c not in sums[n-1]:
    print('no solutions')

  def solutions(i, x, remainder):
    if i==0:
      x_ones = [k for k, v in summands[0].items() if v == remainder]
      if len(x_ones)!=1:
        print('ERROR (internal): wrong sum')
        print(summands[0], remainder, x_ones)
      x[0] = x_ones[0]
      yield x
    else:
      for xi, si in summands[i].items(): # si = ai*xi*xi
        if remainder - si in sums[i-1]:
          x[i]=xi
          yield from solutions(i-1, x.copy(), remainder - si)
  yield from solutions(n-1,[0]*n,c)
  

if __name__ == '__main__':
  q = [1,2,3,4]
  offset = [0,2,3,4]
  print([s for s in squares_sum_solve(q, 2)])
  print([s for s in squares_sum_solve(q, 2, offset = offset)])