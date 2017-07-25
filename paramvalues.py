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
                raise Exception("Value of variable not set.")
                # self.params[param] = random.random() / 1000
            values[param] = self.params[param]
        # print values
        return values

    def set_params_values(self, params):
        for param in params:
            self.params[param] = params[param]

                