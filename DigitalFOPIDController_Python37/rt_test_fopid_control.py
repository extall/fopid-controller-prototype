from fopid_control import fpid
from struct import pack, unpack
import socket

'''
This is a proof-of-concept code for running real-time
FO PID control from Python using a UDP socket
'''

UDP_IP = "127.0.0.1"
UDP_PORT_REMOTE = 5101
UDP_PORT_LOCAL = 5110
UDP_BUFFER_SIZE = 2048

# Set the FOPID parameters. The parameters of the approximation are NOT set in this example
params = {"fpid": {"Kp": -0.002934, "Ki": 0.01030, "Kd": 0.05335, "lam": 0.9, "mu": 0.5}}

# This initializes and computes the controller. From now on, you can access all of the data
# of this approximation *and* can run the control algorithm.
fpid_c = fpid.fpid_2iir(params)

print("Starting control system server...")
print ("UDP local IP:", UDP_IP)
print ("UDP local port:", UDP_PORT_LOCAL)
print ("UDP remote port:", UDP_PORT_REMOTE)

locsock = socket.socket(socket.AF_INET, # Internet
                        socket.SOCK_DGRAM) # UDP

remsock = socket.socket(socket.AF_INET, # Internet
                        socket.SOCK_DGRAM) # UDP

locsock.bind((UDP_IP, UDP_PORT_LOCAL))

# Run the control loop
while True:
    data, addr = locsock.recvfrom(UDP_BUFFER_SIZE)

    # Decode data
    inp = unpack("<d", data)[0]
	
    # Run the control algorithm
    out = fpid_c.compute_fopid_control_law(inp)

    # Send the control law
    msg = pack("<d", out)
    remsock.sendto(msg, (UDP_IP, UDP_PORT_REMOTE))
	
    print("Recieved:{0}, {1} \t\t sent:{2}, {3}".format(inp, data, out, msg))