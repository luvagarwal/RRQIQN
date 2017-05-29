from collections import defaultdict
import random


class Singleton(type):
    _instance = None
    def __call__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instance

class ParamsSingleton(object):
    __metaclass__ = Singleton

    def __init__(self):
        self.params = defaultdict(lambda: None)

    def get_params_values(self, params):
        values = {}
        for param in params:
            if self.params[param] is None:
                self.params[param] = random.random()
            values[param] = self.params[param]

        return values
