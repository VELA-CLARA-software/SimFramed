import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../../../OCELOT'))
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from ocelot.adaptors.genesis import run_genesis
import zmq
import threading
import signal

class genesis_client(object):

    def __init__(self, port=9875, host='apclara1.dl.ac.uk'):
        super(genesis_client, self).__init__()
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REQ)
        self.socket.linger = 0
        self.socket.connect ("tcp://%s:%s" % (host, port))

    def run_genesis(self, inp, launcher, *args, **kwargs):
        self.socket.send_pyobj(inp, launcher, args, kwargs)
        g = self.socket.recv_pyobj()
        return g
