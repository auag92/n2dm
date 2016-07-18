"""A subject base class"""

class Subject:
    """A subject base class"""
    def __init__(self):
        self.observers = []
        self.nsteps = 0

    def attach(self, function, interval=1, *args, **kwargs):
        """Attach callback function.

        At every *interval* steps, call *function* with arguments
        *args* and keyword arguments *kwargs*."""

        if not callable(function):
            function = function.write
        self.observers.append((function, interval, args, kwargs))

    def call_observers(self):
        self.nsteps += 1
        for function, interval, args, kwargs in self.observers:
            if self.nsteps % interval == 0:
                function(*args, **kwargs)

