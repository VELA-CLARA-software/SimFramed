import os, sys
import time
import zmq
import threading
import signal
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.ion()
plt.show()
# fig = plt.figure()

def plotting(alldata, bestdata):
    plt.clf()
    # plt.ion()
    newlist = list(sorted(alldata,key=lambda l: l[0]))
    x,y = zip(*newlist)
    plt.plot(x,y)
    newlist = list(sorted(bestdata,key=lambda l: l[0]))
    x,y = zip(*newlist)
    plt.plot(x,y)
    plt.pause(0.001)
    return None

class zmqServer(threading.Thread):

    _id = -1

    def __init__(self, port=5558):
        super(zmqServer, self).__init__()
        self.event = threading.Event()
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REP)
        self.socket.linger = 0
        self.port = port
        self._best = 1e12
        self._alldata = []
        self._allbest = []

    def run(self):
        self.socket.bind("tcp://*:%s" % (self.port))
        self.socket.linger = 0
        while not self.event.is_set():
            [msg, value, id] = self.socket.recv_pyobj()
            # print ('msg = ', msg)
            if msg == 'get_best':
                self.socket.send_pyobj(self._best)
            elif msg == 'reset_best':
                self._alldata = []
                self._allbest = []
                self._best = 1e12
                self.socket.send_pyobj(self._best)
            elif msg == 'set_best':
                # print('value = ', value)
                if value < self._best:
                    self._best = value
                    print('new best! ', self._best)
                self._allbest.append([id, self._best])
                self._alldata.append([id, value])
                # plotting(self._alldata, self._allbest)
                self.socket.send_pyobj(self._best)

    def stop(self):
        print ('stopping server...')
        self.event.set()
        self.socket.close()
        self.context.destroy()

if __name__ == "__main__":
    server = zmqServer()
    server.daemon = True
    server.start()
    while True:
        time.sleep(1)
