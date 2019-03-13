import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../../../OCELOT'))
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from ocelot.adaptors.genesis import run_genesis
import zmq
import threading
import signal

class genesis_server(threading.Thread):

    def __init__(self, port=9875):
        super(zmqServer, self).__init__()
        self.event = threading.Event()
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REP)
        self.socket.linger = 0
        self.port = port

    def run(self):
        self.socket.bind("tcp://*:%s" % (self.port))
        self.socket.linger = 0
        while not self.event.is_set():
            inp, launcher, args, kwargs = self.socket.recv_pyobj()
            g = self.run_genesis()
            self.socket.send_pyobj(g)

    def stop(self):
        print 'stopping server...'
        self.event.set()
        self.socket.close()
        self.context.destroy()

    def run_genesis(self, inp, launcher, *args, **kwargs):
        g = run_genesis(inp,launcher,i_aft=i_aft)
        return g
