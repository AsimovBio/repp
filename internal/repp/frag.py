from enum import Enum
from math import ceil

from config import conf

class FragType(Enum):
    LINEAR = 0
    CIRCULAR = 1
    PCR = 2
    SYNTHETIC = 3

class Primer:
    def __init__(
            self,
            seq,
            strand,
            penalty,
            pair_penalty,
            tm,
            gc,
            range):
        self.seq = seq
        self.strand = strand
        self.penalty = penalty
        self.pair_penalty = pair_penalty
        self.tm = tm
        self.gc = gc
        self.range = range

class Frag:
    def __init__(
            self,
            id,
            seq,
            start,
            end,
            frag_type,
            unique_id=None,
            type=None,
            cost=None,
            url=None,
            pcr_seq=None,
            primers=None,
            full_seq=None,
            db=None,
            feature_start=None,
            feature_end=None,
            assemblies=None):
        self.id = id
        self.type = type
        self.cost = cost
        self.url = url
        self.seq = seq
        self.pcr_seq = pcr_seq
        self.primers = primers
        self.frag_type = frag_type
        self.unique_id = unique_id
        self.full_seq = full_seq
        self.db = db
        self.start = start
        self.end = end
        self.feature_start = feature_start
        self.feature_end = feature_end
        self.assemblies = assemblies
        self.conf = conf

    def dist_to(self, other):
        return other.start - self.end
    
    def overlaps_via_pcr(self, other):
        return self.dist_to(other) <= conf.pcr_primer_max_embed_length
    
    def overlaps_via_homology(self, other):
        return self.dist_to(other) <= conf.fragments_min_junction_length
    
    def synth_dist(self, other):
        if self.overlaps_via_pcr(other):
            return 0

        clamped_dist = max(1, self.dist_to(other))
        return int(ceil(clamped_dist / conf.synthetic_max_length))

    def cost_to(self, other):
        needs_pcr = self.frag_type in (FragType.PCR, FragType.CIRCULAR)
        pcr_no_homology = 50 * conf.pcr_bp_cost
        pcr_homology = (50 + conf.fragments_min_junction_length) * conf.pcr_bp_cost

        if self == other:
            if needs_pcr:
                return pcr_no_homology
            else:
                return 0
        
        if self.overlaps_via_pcr(other):
            if self.overlaps_via_homology(other):
                return pcr_no_homology
            else:
                return pcr_homology
        
        dist = self.dist_to(other) + conf.fragments_min_junction_length * 2
        frag_cost = synth_fragment_cost(dist)

        if needs_pcr:
            return frag_cost * pcr_no_homology
        else:
            return frag_cost

    def reach(self, nodes, i, features):
        reachable = []

        while True:
            i += 1

            if i > len(nodes):
                return reachable
            if features and nodes[i].feature_end <= self.feature_end:
                continue
            elif nodes[i].end < self.end:
                continue
            
            reachable.append(i)
            
            if nodes[i].unique_id == self.unique_id:
                break

        return

    def junction(self, other, min_homology, max_homology):
        seq1 = (self.pcr_seq or self.seq).upper()
        seq2 = (other.pcr_seq or other.seq).upper()
        start = min(0, len(seq1) - max_homology)
        end = min(0, len(seq2) - min_homology)

        for i in range(start, end):
            sub_seq = seq1[i:]
            if seq2.startswith(sub_seq):
                return sub_seq

            # k = i
            # j = 0
            # while k < len(seq1) and j < len(seq2) and seq1[k] == seq2[j]:
            #     if k == len(seq1) - 1:
            #         return seq1[i:]
            #     j += 1
            #     k += 1

    def synth_to(other, target):
        junction_length = conf.fragments_min_junction_length

        frag_count = self.synth_dist(other)
        if frag_count == 0:
            return None

        total_length = len(target)
        frag_length = max(
            conf.synthetic_min_length,
            (self.dist_to(other) / frag_count) + (junction_length * 2))

        # guessing this is to handle selecting from a plasmid, so we want to wrap
        # around to the beginning again? not sure why it's 4 and not 2 though
        target = (target * 4).upper()
        synths = []
        start = self.end - junction_length + total_length

        while len(synths) < frag_count:
            end = start + frag_length + 1
            seq = target[start:end]

            # TODO: check for hairpins

            synths.append(
                Frag(
                    f'{self.id}-{other.id}-synthesis-{len(synths)+1}',
                    seq,
                    start,
                    end,
                    FragType.SYNTHETIC))

            start = end - junction_length
    
        return synths

    