r"""
Layout for flat flipper triangulations in Sage.
"""

from sage.structure.sequence import Sequence

from sage.categories.fields import Fields

from sage.rings.all import ZZ, QQ, AA, RDF, NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector

from sage.misc.prandom import shuffle
from sage.plot.graphics import Graphics
from sage.plot.line import line2d
from sage.plot.polygon import polygon2d
from sage.plot.text import text
from sage.plot.bezier_path import bezier_path

from .constants import BLUE, RED, HORIZONTAL, VERTICAL
from .misc import flipper_nf_to_sage, flipper_nf_element_to_sage
from .triangulation import Triangulation

QQx = PolynomialRing(QQ, 'x')
_Fields = Fields()

def edge_label(i):
    if i< 0:
        return '~%d' % (-i-1)
    else:
        return str(i)
    
def draw_my_segment(s, t, pos, **kwds):
    x,y = (pos[t] - pos[s])
    if y.is_zero():
        col = 'purple'
    elif x.is_zero():
        col = 'green'
    elif x.sign() == y.sign():
        col = 'red'
    else:
        col = 'blue'
    return line2d([pos[s],pos[t]], color=col, **kwds)

def draw_my_triangle(e1,e2,e3,pos):
    G = Graphics()

    G += polygon2d([pos[e1],pos[e2],pos[e3],pos[e1]], alpha=0.41, color='gray')
    
    G += draw_my_segment(e1, e2, pos)
    G += draw_my_segment(e2, e3, pos)
    G += draw_my_segment(e3, e1, pos)

    G += text(edge_label(e1), (8*pos[e1]+5*pos[e2]+pos[e3])/14, color='black')
    G += text(edge_label(e2), (8*pos[e2]+5*pos[e3]+pos[e1])/14, color='black')
    G += text(edge_label(e3), (8*pos[e3]+5*pos[e1]+pos[e2])/14, color='black')

    return G

def orientation(u, v, w):
    a = v - u
    b = w - u
    return (a[0]*b[1] - a[1]*b[0]).sign()

def has_intersection(triangles, new_triangle, pos):
    r"""
    Check whether ``new_triangle`` intersects one of the triangles in ``triangles``
    where the positions of the vertices have to be found in ``pos``.
    """
    for t in triangles:
        for i in range(3):
            p1 = t[i]
            p2 = t[(i+1)%3]
            for j in range(3):
                q1 = new_triangle[j]
                q2 = new_triangle[(j+1)%3]

                if (orientation(p1,p2,q1) != orientation(p1,p2,q2) and \
                    orientation(q1,q2,p1) != orientation(q1,q2,p2)):
                    return True

    return False

class FlatTriangulation:
    def __init__(self, triangles, vectors):
        r"""
        INPUT: triangle (or triangulation) and vectors
        """
        if isinstance(triangles, Triangulation):
            triangles = triangles.faces()

        self._triangulation = Triangulation(triangles)
        n = self._triangulation.num_edges()
        self._edge_to_face = [None] * (2*n)
        self._faces = self._triangulation.faces()

        for i,(e0,e1,e2) in enumerate(self._triangulation):
            self._edge_to_face[e0] = i
            self._edge_to_face[e1] = i
            self._edge_to_face[e2] = i

        if len(vectors) == n:
            vectors = list(vectors)
            vectors.extend(vectors[::-1])
        if len(vectors) != 2*n:
            raise ValueError('wrong number of vectors')
        
        vectors = Sequence([vector(v) for v in vectors])
        self._V = vectors.universe()
        self._K = self._V.base_ring()

        if self._K not in _Fields:
            self._K = self._K.fraction_field()
            self._V = self._V.change_ring(self._K)
            vectors = [v.change_ring(self._K) for v in vectors]

        self._vectors = list(vectors)

        # small check
        for a in range(n):
            u = self._vectors[a]
            v = self._vectors[~a]
            if u != v and u != -v:
                raise ValueError('vec[%s] = %s while vec[%s] = %s' % (edge_label(a), u, edge_label(~a), v))

        for a,b,c in self._triangulation:
            va = self._vectors[a]
            vb = self._vectors[b]
            vc = self._vectors[c]
            if va + vb + vc:
                raise ValueError('vec[%s] = %s, vec[%s] = %s and vec[%s] = %s do not sum to zero' % (edge_label(a), va, edge_label(b), vb, edge_label(c), vc))

        # spanning tree for display
        self._root = None
        self._edges = None

    @staticmethod
    def from_flipper(F):
        triangles = [(e.label, f.label, g.label) for e,f,g in F.triangulation]
        n = 3 * len(triangles) / 2
        x = F.edge_vectors.values()[0][0]
        K = flipper_nf_to_sage(x.number_field)
        V = VectorSpace(K, 2)
        # translate into Sage number field
        Vec = {e.label: V((flipper_nf_element_to_sage(v.x, K),
                               flipper_nf_element_to_sage(v.y, K)))
                  for e,v in F.edge_vectors.items()}
        vectors = [Vec[i] for i in range(n)]
        vectors.extend(Vec[i] for i in range(-n, 0))

        return Triangulation(triangles, vectors)

    @staticmethod
    def from_coloured_triangulation(T):
        PH = T.train_track_polytope(HORIZONTAL)
        PV = T.train_track_polytope(VERTICAL)

        # normalization
        PH.add_constraint(c)
        PV.add_constraint(c)

        # pick vertices
        PH.generators()
        PV.generators()

    def xy_scaling(self, sx, sy):
        for i,v in enumerate(self._vectors):
            self._vectors[i] = self._V((sx*v[0], sy*v[1]))

    def __repr__(self):
        return 'Flat Triangulation made of %s triangles' % len(self._triangulation)

    def set_random_spanning_tree(self, start=0):
        # just select the edges
        start = int(start)
        self._root = start
        n = len(self._triangulation)
        face_seen = [False] * n
        faces = []

        wait = [start]
        face_seen[start] = True
        faces.append(start)

        edges = self._edges = {t:[] for t in range(n)}
        fp = self._triangulation.face_permutation()
        while wait:
            shuffle(wait)
            t = wait.pop()
            # choice of edge is not random here!!!
            e0 = fp[t]
            t_edges = list(self._faces[t])
            shuffle(t_edges)
            for e in t_edges:
                assert self._edge_to_face[e] == t
                f = self._edge_to_face[~e]
                if face_seen[f] is False:
                    edges[t].append(e)
                    face_seen[f] = True
                    faces.append(f)
                    wait.append(f)
        
    def _plot_train_track(self, slope, pos, For, color):
        V2 = VectorSpace(RDF, 2)
        G = Graphics()

        if slope == HORIZONTAL:
            POS = RED
            NEG = BLUE
        else:
            POS = BLUE
            NEG = RED

        for e0, e1, e2 in self._triangulation:
            # determine the large edge
            v0 = For[e0]
            v1 = For[e1]
            v2 = For[e2]
            col0 = RED if v0[0] * v0[1] > 0 else BLUE
            col1 = RED if v1[0] * v1[1] > 0 else BLUE
            col2 = RED if v2[0] * v2[1] > 0 else BLUE
            if col0 == POS and col1 == NEG:
                # e2 is large
                l,s1,s2 = e2,e0,e1
            elif col1 == POS and col2 == NEG:
                # i is large
                l,s1,s2 = e0,e1,e2
            elif col2 == POS and col0 == NEG:
                # j is large
                l,s1,s2 = e1,e2,e0

            pl = V2(pos[l])
            ps1 = V2(pos[s1])
            ps2 = V2(pos[s2])

            cl = (pl + ps1) / 2
            vl = (ps1 - pl)
            cs1 = (ps1 + ps2) / 2
            vs1 = ps2 - ps1
            cs2 = (ps2 + pl) / 2
            vs2 = pl - ps2

            ol = V2((-vl[1], vl[0]))
            ol /= ol.norm()
            os1 = V2((-vs1[1], vs1[0]))
            os1 /= os1.norm()
            os2 = V2((-vs2[1], vs2[0]))
            os2 /= os2.norm()

            G += bezier_path([[cl, cl + 0.3 * ol,
                               cs1 + 0.3 * os1, cs1]], rgbcolor=color)
            G += bezier_path([[cl, cl + 0.3 * ol,
                               cs2 + 0.3 * os2, cs2]], rgbcolor=color)
        return G


    def plot(self, horizontal_train_track=False,
                   vertical_train_track=False):
        if self._root is None:
            self.set_random_spanning_tree()

        n = len(self._triangulation)

        For = self._vectors[:]
        pos = [None] * (3*n)   # edge (actually half-edge) position is much clever !
        G = Graphics()

        # set the position of the first triangle
        e1,e2,e3 = self._faces[self._root]
        pos[e1] = self._V.zero()
        pos[e2] = pos[e1] + For[e1]
        pos[e3] = pos[e2] + For[e2]
        G += draw_my_triangle(e1,e2,e3,pos)

        wait = self._edges[self._root][:] # edge waiting set

        while wait:
            e = wait.pop()
            t = self._edge_to_face[~e]  # triangle number seen trough e
            f = self._faces[t]          # triangle seen through e

            # compute edges in the right order
            e1 = ~e
            i = f.index(e1)
            e2 = f[(i+1)%3]
            e3 = f[(i+2)%3]

            # possibly fix orientation
            if For[~e] != -For[e]:
                For[e1] *= -1
                For[e2] *= -1
                For[e3] *= -1

            pos[e1] = pos[e] + For[e]
            pos[e2] = pos[e1] + For[e1]
            pos[e3] = pos[e2] + For[e2]
            G += draw_my_triangle(e1,e2,e3,pos)
            wait.extend(self._edges[t])

        if horizontal_train_track:
            G += self._plot_train_track(HORIZONTAL, pos, For, 'green')
        if vertical_train_track:
            G += self._plot_train_track(VERTICAL, pos, For, 'purple')

        #print pos
        G.set_aspect_ratio(1)
        G.axes(False)
        return G
