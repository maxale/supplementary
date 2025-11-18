'''
Supplementary SageMath code ver. 20251118 for the paper:
* N. Alexeev and M. A. Alekseyev. "Combinatorial Scoring of Phylogenetic Trees and Networks Based on Homoplasy-Free Characters". Journal of Computational Biology 25:11 (2018), 1203-1219.
  DOI: https://doi.org/10.1089/cmb.2018.0082
  Preprint: arXiv:1602.02841 [q-bio.PE] https://doi.org/10.48550/arXiv.1602.02841
'''

def FGH( N, v, true_leaves = None ):
  '''
  Input:
    `N`, a cactus network, possibly a subnetwork of the bigger ("parent") network we started to work with
    `v`, a vertex in `N`
    `true_leaves`, a set of leaves in the parent network; it is omitted when `N` is a parent network
  Output: 
    [ F_v(x), G_v(x), H_v(x) ]
  '''

  #print( "Input: ",N.edges(),"\t",v )

  if not N.is_directed_acyclic():
    raise ValueError("Input graph is not a network!")

  if true_leaves is None:
    # N is a parent network, initialize true_leaves:
    true_leaves = Set(N.sinks())

  # define g.f. variable x
  x = QQ['x'].0

  if N.out_degree(v) == 0:
    # v is a leaf
    #print( v," is a leaf" )
    if v in true_leaves:
      #print( "Output: ",[x,0,0] )
      return [x,0,0]
    else:
      #print( "Output: ",[0,0,1] )
      return [0,0,1]

  # Notice that posets 'grow' up
  P = Poset(N)

  # looking for a sink t
  t = v
  S = [Set(P.principal_upper_set(u)) for u in N.neighbors_out(v)]
  for u in Combinations(S,2):
    T = u[0].intersection(u[1])
    if not T.is_empty():
      # sink is found
      t = P.subposet(T).minimal_elements()[0]
      break

  if t == v:
    # v is a regular vertex
    #print( v," is regular" )

    # compute FGH at the children of v
    cFGH = [FGH(N,u,true_leaves) for u in N.neighbors_out(v)]

    H = prod(fgh[0]+fgh[2] for fgh in cFGH)
    G = sum((fgh[0]+fgh[1])*H/(fgh[0]+fgh[2]) for fgh in cFGH)
    F = x*prod(fgh[0]+fgh[2] + (fgh[0]+fgh[1])/x for fgh in cFGH) - x*H - G;

  else:
    # v is a source
    #print( v," is a source having a sink ",t )

    # F_t, G_t, H_t
    t_FGH = FGH(N, t, true_leaves)

    p = N.neighbors_in(t)
    if len(p) != 2:
      raise ValueError("Parents error p")

    N.delete_edge(p[0],t)
    # N represents L at this point
    L_FGH = FGH(N, v, true_leaves)
    N.add_edge(p[0],t)

    N.delete_edge(p[1],t)
    # N represents R at this point
    R_FGH = FGH(N, v, true_leaves)

    N.delete_edge(p[0],t)
    # N represents N_{s\t} U N_t at this point
    Nst_FGH = FGH(N, v, true_leaves)

    # remove edges in the branching path v => t from N to form B
    for q in p:
      u = q
      while u!=v:
        w = N.neighbors_in(u)
        if len(w)!=1:
           raise ValueError("Parents error w")
        N.delete_edge(w[0],u)
        u = w[0]

    # N represents B at this point
    B_FGH = [FGH(N, u, true_leaves) for u in Set(P.closed_interval(v,p[0])) + Set(P.closed_interval(v,p[1]))]

    B_prodH = prod(fgh[2] for fgh in B_FGH)

    H = L_FGH[2] + R_FGH[2] - Nst_FGH[2] * (t_FGH[0] + t_FGH[2])

    G = L_FGH[1] + R_FGH[1] - Nst_FGH[1] * (t_FGH[0] + t_FGH[2]) - (t_FGH[0] + t_FGH[1]) * B_prodH

    F = L_FGH[0] + R_FGH[0] - Nst_FGH[0] * (t_FGH[0] + t_FGH[2]) - (t_FGH[0] + t_FGH[1]) * (prod((fgh[0]+fgh[1])/x + fgh[2] for fgh in B_FGH) - B_prodH)

  #print( "Output: ",[F,G,H] )
  return [F,G,H]


def Fig5():

  fgh = FGH( DiGraph( {6:[1,5,7,9], 7:[2,8], 9:[4,8], 8:[3]} ), 6 )
  print( fgh )
  print( [fgh[0](1), fgh[1](1), fgh[2](1)] )

  fgh = FGH( DiGraph( {6:[1,5,7,9], 7:[2,8], 9:[4], 8:[3,9]} ), 6 )
  print( fgh )
  print( [fgh[0](1), fgh[1](1), fgh[2](1)] )
