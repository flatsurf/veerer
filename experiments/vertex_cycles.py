
from itertools import product
from veerer import *
from surface_dynamics import *
from queue import Queue

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
    
    current = vertex_cycles(T1)
    future = set(flip.inverse()(curve) for curve in vertex_cycles(T2))
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
    current, future = current_future_vertex_cycles(T, edge, colour)
    future_twists = set(curve.encode_twist(power=k) for curve in future for k in [-1, 1])
    
    
    to_do = Queue()
    to_do.put(T.to_curver().id_encoding())
    seen = {T.to_curver().id_encoding(): 0}
    
    for cycle in current:
        if cycle in future: continue
        
        twist = cycle.encode_twist()
        while twist not in seen:
            g = to_do.get()
            for h in future_twists:
                f = g * h
                if f not in seen:
                    seen[f] = seen[g] + 1
                    to_do.put(f)
                    print(seen[f], str(f.self_image()))
                    if len(seen) % 100 == 0: print('expanding seen past {}'.format(len(seen)))

def test(T, edge, colour):
    convert = lambda X: (X[0], tuple(X[1].flatten()))  # Since numpy.ndarrays are not hashable we need a converter.
    current, future = current_future_vertex_cycles(T, edge, colour)
    twists = [(curve, k, curve.encode_twist(power=k), curve.encode_twist(power=k).homology_matrix()) for curve in future for k in [-1, +1]]
    
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
            for curve, k, t, M in twists:
                next_lam_mat = (t(lam), M.dot(mat))
                key = convert(next_lam_mat)
                if key not in depths:
                    depths[key] = [(curve, k)] + depths[convert(current_lam_mat)]
                    to_do.put(next_lam_mat)
                    if len(depths) % 10000 == 0: print('expanding seen past {}'.format(len(depths)))
        print(depths[twist_key])



def test_all_in(stratum):
    for T in Automaton.from_stratum(stratum):
        print(T)
        assert T.is_core()
        for edge, colour in core_flips(T):
            print(T, edge, colour)
            test(T, edge, colour)

def intersection_matrix(T):
    ''' Return the matrix of intersection number of vertex cycles of T. '''
    cycles = list(vertex_cycles(T))
    return [[x.intersection(y) for y in cycles] for x in cycles]

if __name__ == '__main__':
    # T = VeeringTriangulation("(0,~2,1)(2,4,~3)(3,~5,~4)(5,~1,~0)", "RBBRBB")
    # test(T, 2, BLUE)
    test_all_in(AbelianStratum(0, 0, 0))

