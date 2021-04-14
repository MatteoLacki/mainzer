import time

class Timer(dict):
    def __init__(self):
        self.time = time.time()

    def __call__(self, name):
        new_time = time.time()
        self[name] = new_time - self.time
        self.time = new_time