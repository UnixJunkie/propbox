###### A simple Future-like implementation

# In the long-term, I think Python 3.4's asyncio is the right
# framework to handle those dependencies. This project is a stepping
# stone to get here. For now I only ue the 'future' API.

PENDING = 'PENDING'
FINISHED = 'FINISHED'

class Future(object):
    def __init__(self):
        self._state = PENDING
        self._exception = None
        self._result = None

    def cancel(self):
        # I only use complete Futures
        return False

    def cancelled(self):
        # Can't cancel these futures
        return False

    def running(self):
        # Public futures are always finished
        return False

    def done(self):
        return True

    def result(self, timeout=None):
        if self._state == FINISHED:
            if self._exception is not None:
                raise self._exception
            else:
                return self._result
            
        raise AssertionError("Nothing computed")

    def exception(self, timeout=None):
        if self._state == FINISHED:
            return self._exception
        raise AssertionError("Nothing computed")

    def add_done_callback(self, fn):
        if self._state == FINISHED:
            fn(self)
            return
        raise AssertionError("Should never see non-done futures")
    
    def set_result(self, result):
        if self._state != PENDING:
            raise AssertionError("Already set; cannot modify")
        self._result = result
        self._state = FINISHED

    def set_exception(self, exception):
        if self._state != PENDING:
            raise AssertionError("Already set; cannot modify")
        self._exception = exception
        self._state = FINISHED
            
# Some helper functions to make a future with a result or exception
        
def new_future(value):
    f = Future()
    f.set_result(value)
    return f

def new_future_exception(exception):
    f = Future()
    f.set_exception(exception)
    return f
