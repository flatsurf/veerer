
from itertools import product
from veerer import *
from surface_dynamics import *
from queue import Queue

import ppl

def compose(P):
    assert P
    h = None
    for p in P:
        h = p * h
    return h

def vertex_cycles(T):
    ''' Return set of vertex cycles as curves on the (flipper) triangulation. '''
    polytope = T.train_track_polytope()
    rays = [gen for gen in polytope.generators() if gen.is_ray()]
    T_flipper = T.to_curver()
    return set(T_flipper.lamination([int(x) for x in ray.coefficients()]) for ray in rays)

def current_future_vertex_cycles(T1, edge, colour):
    assert edge in T1.forward_flippable_edges()
    T2 = T1.copy()
    T2.flip(edge, colour)
    
    T1_flipper = T1.to_curver()
    flip = T1_flipper.encode_flip(~edge)
    
    current = set(flip(curve) for curve in vertex_cycles(T1))
    future = vertex_cycles(T2)
    return current, future

def core_flips(T):
    for edge in T.forward_flippable_edges():
        T0 = T.copy()
        T0.flip(edge, RED)
        if T0.edge_has_curve(edge):
            yield edge, RED
        T0 = T.copy()
        T0.flip(edge, BLUE)
        if T0.edge_has_curve(edge):
            yield edge, BLUE

def test(T, edge, colour):
    ''' Test whether the vertex cycles of T can be written as a product of twists in the vertex cycles of T after the given edge has been split. '''
    convert = lambda X: (X[0], tuple(X[1].flatten()))  # Since numpy.ndarrays are not hashable we need a converter.
    current, future = current_future_vertex_cycles(T, edge, colour)
    assert False
    twists = [(('' if k == 1 else '~') + str(index), curve.encode_twist(power=k), curve.encode_twist(power=k).homology_matrix()) for index, curve in enumerate(future) for k in [-1, +1]]
    
    identity = T.to_curver().id_encoding()
    identity_lam_mat = (identity.self_image(), identity.homology_matrix())
    to_do = Queue()
    to_do.put(identity_lam_mat)
    depths = {convert(identity_lam_mat): []}
    
    for cycle in current:
        twist = cycle.encode_twist()
        twist_key = convert((twist.self_image(), twist.homology_matrix()))
        while twist_key not in depths:
            current_lam_mat = to_do.get()
            lam, mat = current_lam_mat
            for s, t, M in twists:
                next_lam_mat = (t(lam), M.dot(mat))
                key = convert(next_lam_mat)
                if key not in depths:
                    depths[key] = [s] + depths[convert(current_lam_mat)]
                    to_do.put(next_lam_mat)
                    if len(depths) % 10000 == 0: print('expanding seen past {}'.format(len(depths)))
        print('{:s} decomposes as: {}'.format(list(cycle), '.'.join(depths[twist_key])))

def test_conjugators(T, edge, colour):
    ''' Test whether new twists can be applied to the vertex cycles of T to make them carried after the split.
    
    This is the same as seeing whether the twists about the vertex cycles of T can be written as conjugate of a product of twists in the vertex cycles of T'. '''
    current, future = current_future_vertex_cycles(T, edge, colour)
    future = list(future)
    T2 = T.copy()
    T2.flip(edge, colour)
    polytope = T2.train_track_polytope()
    twists = [((index+1) * (+1 if k == 1 else -1), curve.encode_twist(power=k)) for index, curve in enumerate(future) for k in [-1, +1]]
    print('New vertex cycles:')
    for index, curve in enumerate(future):
        print('{}: {}'.format(index+1, list(curve)))
    print('')
    
    for cycle in current:
        if polytope.contains(ppl.polyhedron.C_Polyhedron(ppl.point(ppl.Linear_Expression(list(cycle), 0)))):
            print('{} is already carried'.format(list(cycle)))
            continue
        
        names = {cycle: tuple()}
        to_do = Queue()
        to_do.put(cycle)
        found = False
        while not found:
            current = to_do.get()
            for symbol, twist in twists:
                new = twist(current)
                if new not in names:
                    names[new] = (symbol,) + names[current]
                    to_do.put(new)
                    if len(names) % 10000 == 0: print('expanding seen past {}'.format(len(names)))
                    
                    if polytope.contains(ppl.polyhedron.C_Polyhedron(ppl.point(ppl.Linear_Expression(list(new), 0)))):
                        print('{} became carried after applying {}'.format(list(cycle), '.'.join(map(str, names[new]))))
                        found = True
                        break

def test_all_in(stratum):
    for T in Automaton.from_stratum(stratum):
        print(T)
        assert T.is_core()
        for edge, colour in core_flips(T):
            print(T, edge, colour)
            test_conjugators(T, edge, colour)

def intersection_matrix(T):
    ''' Return the matrix of intersection number of vertex cycles of T. '''
    cycles = list(vertex_cycles(T))
    return [[x.intersection(y) for y in cycles] for x in cycles]

if __name__ == '__main__':
    # T = VeeringTriangulation("(0,~2,1)(2,4,~3)(3,~5,~4)(5,~1,~0)", "RBBRBB")
    # test(T, 2, BLUE)
    # test_all_in(QuadraticStratum(1, 1, -1, -1))
    test_all_in(QuadraticStratum(1, 1, -1, -1))
    # test_all_in(AbelianStratum(0, 0, 0))

