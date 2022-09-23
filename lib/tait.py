import enum
import itertools


class Elevation(enum.Enum):
    UNDER = -1
    OVER = 1


class Side(enum.Enum):
    LEFT = -1
    RIGHT = 1


class Alignment(enum.Enum):
    AGAINST = -1
    WITH = 1


class Sign(enum.ENUM):
    NEG = -1
    POS = 1


class CircleSide(enum.Enum):
    INSIDE = -1
    OUTSIDE = 1


class CrossingItem():

    def __init__(self, crossing: int, elevation: Elevation, sign: Sign):
        assert isinstance(crossing, int)
        assert isinstance(elevation, Elevation)
        self.crossing = crossing
        self.elevation = elevation
        self.sign = sign


class Strand():

    def __init__(self, crossings: tuple[CrossingItem], alignment):
        self.crossings = crossings
        self.alignment = alignment

    def __eq__(self, other):
        if self.crossings == other.crossings and\
                self.alignment == other.alignment:
            return True
        elif self.crossings == reversed(other.crossings) and\
                self.alignment == -other.alignment:
            return True
        else:
            return False


class FaceEdge():

    def __init__(self, strand: Strand, side: CircleSide):
        self.strand = strand
        self.side = side


def tait_graphs(gauss_code: tuple[Strand]):

    def next_face_edge(strand, alignment, turn, jump):
        

    for turn in (Side.LEFT, Side.RIGHT):
        counter = itertools.count()
        crossing_edges = {i: [None, None] for i in len(gauss_code)}
        face_edge_to_face = {}

        discovered = set()
        discovered.add(gauss_code[0])
        while discovered:
            d = discovered.pop()
            tait_vertex = next(counter)
            
