import time
import zmq
import json

class zmqClient(object):

    def __init__(self, port=5558, host='localhost'):
        super(zmqClient, self).__init__()
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REQ)
        self.socket.linger = 0
        self.socket.connect ("tcp://%s:%s" % (host, port))

    def request(self, msg='get_best', value=None, id=None):
        self.socket.send_pyobj([msg, value, id])
        message = self.socket.recv_pyobj()
        # print "Received reply ", kwargs, "[", message, "]"
        return message

    def get_best(self):
        return self.request(msg='get_best')

    def reset_best(self):
        return self.request(msg='reset_best')

    def set_best(self, best, id=None):
        return self.request(msg='set_best', value=best, id=id)

if __name__ == "__main__":
    client = zmqClient()
    print (client.set_best(12,1))
    print (client.set_best(11,2))
    print (client.set_best(10,3))
    print (client.set_best(17,4))
    exit()
