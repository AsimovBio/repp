from os.path import expanduser
from yaml import safe_load
from math import ceil, floor

class SynthCost:
    def __init__(self, fixed, cost):
        self.fixed = fixed
        self.cost = cost

class Config:
    def __init__(self):
        self._config = {}

    def load(self, config_path):
        with open(config_path) as f:
            self._config = safe_load(f)

    def __getattr__(self, attr):
        # convert the attribute to kebab case and look it up in the config
        key = attr.replace('_', '-')
        if key not in self._config:
            raise KeyError(attr)

        return self._config[key]

    def synth_fragment_cost(self, frag_length):
        frag_count = ceil(frag_length / self.synthetic_max_length)
        frag_length = floor(frag_length / frag_count)
        cost = self.synth_cost(frag_length, self.synthetic_fragment_cost)
        if cost.fixed:
            return frag_count * cost.cost
        else:
            return frag_count * frag_length * cost.cost

    def synth_plasmid_cost(self, insert_length):
        cost = self.synth_cost(insert_length, self.synthetic_plasmid_cost)
        if cost.fixed:
            return cost.cost
        else:
            return insert_length * cost.cost

    def synth_cost(self, seq_length, costs):
        # python dictionaries are guaranteed to be ordered (3.6+)
        for length in costs:
            if length >= seq_length:
                # wrap it in a SynthCost object
                return SynthCost(**costs[length])

        return SynthCost(True, float('inf'))

# this is a singleton, change it here
conf = Config()
conf.load(expanduser('~/.repp/config.yaml'))
