#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
from subprocess import call

def gram_matrix2graph(M, d):
    def weight(M, i, j):
        cos2 = (M[i][j]*M[i][j])/(M[i][i]*M[j][j])
        if cos2 == 0:
            return 2
        if cos2 == 1/4:
            return 3
        if cos2 == 1/2:
            return 4
        if cos2 == 3/4:
            return 6
        if cos2 == 1:
            return 0
        if cos2 > 1:
            return 1
        else:
            print('coxiter.py ERROR: cosine ', cos2)
            return -1

    n = len(M)
    strings = ['{0} {1}\n'.format(n,d-1),]
    for i in range(n):
        for j in range(i):
            if M[i][j] != 0:
                strings.append('{0} {1} {2}\n'.format(j+1, i+1, weight(M,i,j) ))
    strings.append('\n')
    print(strings)
    return strings
    

def run(M, d, graph_file = 'graph.txt', answer_file = 'answer.txt', coxiter_dir = './CoxIter/' ):
    with open(coxiter_dir+graph_file, 'w') as graph:
        graph.writelines(gram_matrix2graph(M, d))
    call('bash -c "{0}coxiter -fv < {0}{1} > {0}{2}"'.format(coxiter_dir, graph_file, answer_file), shell=True)
    with open(coxiter_dir+answer_file, 'r') as out:
        response = out.readlines()
    question = 'Finite covolume'
    answer = [('yes' in s) for s in response if question in s]
    if len(answer) == 0:
        print('"{}" not found, check answer.txt, you may want to rebuild CoxIter'.format(question))
        raise Exception('coxiter.py', 'did not find answer')
    return answer[0]
