class map(dict):
    def __init__(self):
        dict.__init__(self)
        self.ordered_keys = []
    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        self.ordered_keys.append(key)
