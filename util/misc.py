from collections import defaultdict

class keydefaultdict(defaultdict):
    def __missing__(self, key):
        if self.default_factory is not None:
            ret = self[key] = self.default_factory(key)
            return ret
        else:
            raise KeyError(key)
