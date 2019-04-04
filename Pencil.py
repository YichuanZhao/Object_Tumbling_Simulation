
import math

import numpy as np

m = 1.0

r = 0.5

h1 = 3.0

h2 = 0.5

m1 = m*18/19

m2 = m/19

h = (6*(h1)**2 + 12*h1*h2 + 3*(h2)**2)/(12*h1 + 4*h2)

l = h1/2 + h2 - h

Q11 = (m/(h1 + h2/3))*( h1*((3*r**2 + h1**2)/12 + l**2) + (h2/3)*((3/5)*(r**2/4 + h2**2) + h**2))

Q22 = Q11

Q33 = (0.5*m1 + 0.3*m2)*r**2

Q = np.matrix([[Q11, 0, 0],[0, Q22, 0],[0, 0, Q33]])

Q_inverse = np.matrix([[1/Q11, 0, 0],[0, 1/Q22, 0],[0, 0, 1/Q33]])

t1 = 0.4

t2 = 0.8

p1 = h*np.array([math.cos(math.pi/3), 0.0, math.sin(math.pi/3)])

q1 = np.array([math.cos(math.pi/12), 0, math.sin(math.pi/12), 0]) 

v1_n = -5*np.array([math.cos(math.pi/6), 0.0, math.sin(math.pi/6)])

w1_n = np.array([1.0, 5.0, 0.5])

v1_p = np.array([-1.80954, -0.546988, 1.2076])

w1_p = np.array([0.09957, -0.04174, 0.5])

g = np.array([0, 0, -9.8])

# print(m1)
# print(m2)
# print(h)
# print(l)
# print(Q11)
# print(Q22)
# print(Q33)
# print(p1)



def dotProduct(a, b):
	res = 0

	for i in range(len(a)):
		res += a[i] * b[i]	

	return res


test1 = np.matrix([[1,2,3], [1,2,3], [1,2,3]])
test2 = np.array([1,1,1])



def crossProduct(a, b):
	res = np.zeros(3)

	res[0] = a[1]*b[2] - a[2]*b[1]

	res[1] = a[2]*b[0] - a[0]*b[2]

	res[2] = a[0]*b[1] - a[1]*b[0]

	return res

def getVector(M, V):
	res = np.zeros(3)
	for i in range(3):
		res[i] = 0
		for j in range(3):
			res[i] += M[i,j]*V[j]

	return res

def quaProduct(a, b):
	res = np.zeros(4)
	res[0] = a[0] + b[0] - dotProduct(a[1:4], b[1:4])
	res[1:4] = a[0]*b[1:4] + b[0]*a[1:4] + crossProduct(a[1:4], b[1:4])

	return res	

def tumble(p_s, q_s, v_s, w_s, T_s, T_e):

	t = 0

	v = v_s

	w = w_s

	q = q_s

	p = p_s

	T = T_s - T_e

	p_t = 0.00001

	count = 0

	res_p = []

	res_v = []

	res_w = []

	while math.fabs(t) <= math.fabs(T):
		omega = np.zeros(3)
	    
		omega = (q[0]**2 - dotProduct(q[1:4], q[1:4]))*v + 2*dotProduct(q[1:4], v)*q[1:4] + 2*q1[0]*crossProduct(q[1:4], v)

		omega_abs = dotProduct(omega, omega)**(0.5)
		omega_u = omega/omega_abs

		phi_u = 0.00

		if T > 0:
			phi_u = omega_abs*p_t
		else:
			phi_u = -omega_abs*p_t

		R = np.zeros(4)

		R[0] = math.cos(phi_u/2)
		R[1:4] = omega_u*math.sin(phi_u/2)

		q = quaProduct(R, q)

		if T > 0:
			v = v + p_t*g
			w = w - p_t*(getVector(Q_inverse,(crossProduct(w, getVector(Q, w))))) 
			t = t + p_t

		else:
			v = v - p_t*g
			w = w + p_t*(getVector(Q_inverse,(crossProduct(w, getVector(Q, w))))) 
			t = t - p_t

		p = p_s + t*v_s + 0.5*(t**2)*g

		count += 1

		if count%10000 == 0:
			# print('At time {}, the pencil position is {}, velocity is {}, and angular velocity is {}'.format(t, p, v, w))
			res_p.append(p)
			res_v.append(v)
			res_w.append(w)

	return res_p, res_v, res_w


# tumble(p1, q1, v1_p, w1_p, 0.4)

temp_p, temp_v, temp_w = tumble(p1, q1, v1_n, w1_n, 0, 0.4)

for i in range(4):
	print('At time {}, the pencil position is {}, velocity is {}, and angular velocity is {} \n'.format(i*0.1, temp_p[3-i], temp_v[3-i], temp_w[3-i]))


print('At time {}-, the pencil position is {}, velocity is {}, and angular velocity is {}\n'.format(0.4, p1, v1_n, w1_n))

print('At time {}+, the pencil position is {}, velocity is {}, and angular velocity is {} \n'.format(0.4, p1, v1_p, w1_p))


temp_p, temp_v, temp_w = tumble(p1, q1, v1_p, w1_p, 0.8, 0.4)

for i in range(4):
	print('At time {}, the pencil position is {}, velocity is {}, and angular velocity is {} \n'.format((i+5)*0.1, temp_p[i], temp_v[i], temp_w[i]))	


#print(p1 + 0.4*v1_p + 0.5*(0.16)*g)