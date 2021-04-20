# This is just a quick example showing how to configure the control server with new FOPID parameters
import socket
from struct import pack, unpack

UDP_IP = "127.0.0.1"
UDP_PORT_CTRL = 5201

confs = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

# Set the parameters here to test and send them using the socket
Kp = -0.01
Ki = 0.25
Kd = 0.05
lam = 0.8
mu = 0.6

# NB! Parameter order is important! "c" means change FOPID parameters
msg = pack("<cddddd", "c".encode('ascii'), Kp, Ki, Kd, lam, mu)
confs.sendto(msg, (UDP_IP, UDP_PORT_CTRL))