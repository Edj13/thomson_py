#!/bin/python

import sys
import numpy as np
from scipy import optimize
import random
import itertools

# Global Vars
rad = float(1)
radsq = rad*rad
npts = 12

# Functions
def printRes(x,q):
#    print "Results"
    xx = np.column_stack((x[:npts], x[npts:], q[:npts], q[npts:]))
    for atom in xx:
	st1 = np.sin(atom[0])
	ct1 = np.cos(atom[0])
	sp1 = np.sin(atom[1])
	cp1 = np.cos(atom[1])
	x1 = ct1*sp1
	y1 = st1*sp1
	z1 = cp1
	print x1,y1,z1
	
def writeLammps(x,q, writeName):
    g=open(writeName, "a")
    g.write("ITEM:TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n")
    g.write(str(npts)+"\nITEM: BOX BOUNDS pp pp pp\n")
    g.write("-2 2\n-2 2\n-2 2\n")
    g.write("ITEM: ATOMS id type x y z q\n")
    xx = np.column_stack((x[:npts], x[npts:], q[:npts], q[npts:]))
    ctr = 1
    for atom in xx:
	st1 = np.sin(atom[0])
	ct1 = np.cos(atom[0])
	sp1 = np.sin(atom[1])
	cp1 = np.cos(atom[1])
	x1 = ct1*sp1
	y1 = st1*sp1
	z1 = cp1
	g.write(str(atom[3])+" 1 "+str(x1)+" "+str(y1)+" "+str(z1)+" "+str(atom[2])+"\n")

    g.close()





def genPts(n):
    a = np.random.rand(n*2,1)
    a[:npts] = a[:npts]*2*np.pi
    a[npts:] = a[npts:]*np.pi
#    charges = np.zeros(npts*2).reshape(npts,2)
#    charges.fill(1)
#    charges[:,1] = np.arange(0, npts)
    charges = np.zeros(npts*2)
    ones = np.empty(npts)
    ones.fill(1)
    charges[:npts] = ones
    charges[npts:] = np.arange(0, npts)
#    print charges

    return a,charges

def coulomb(x,q):
    fx = 0
    dfx = np.zeros(npts*2)
#    print "x: ",x
#    print "q: ",q
    xx = np.column_stack((x[:npts], x[npts:], q[:npts], q[npts:]))
    for pair in itertools.combinations(xx,2):
	st1 = np.sin(pair[0][0])
	ct1 = np.cos(pair[0][0])
	sp1 = np.sin(pair[0][1])
	cp1 = np.cos(pair[0][1])
	st2 = np.sin(pair[1][0])
	ct2 = np.cos(pair[1][0])
	sp2 = np.sin(pair[1][1])
	cp2 = np.cos(pair[1][1])
	c12 = np.cos(pair[0][1]-pair[1][1])
	s12 = np.sin(pair[0][1]-pair[1][1])
	x1 = ct1*sp1
	y1 = st1*sp1
	z1 = cp1
	x2 = ct2*sp2
	y2 = st2*sp2
	z2 = cp2
#	print "1: ", x1,y1,z1
#	print "2: ", x2,y2,z2
	

	drsq = 2*radsq - 2*radsq*(st1*st2*(cp1*cp2+sp1*sp2)+ct1*ct2)
#	print "st1*st2: ", st1*st2
#	print "ct1*ct2: ", ct1*ct2
#	print "c12: ", c12
#	print "cp1*cp2+sp1*sp2: ", (cp1*cp2+sp1*sp2)
#	print "hi: ",(st1*st2*(cp1*cp2+sp1*sp2)+ ct1*ct2)
#	print "drsq: ", drsq
	drsq = 2*radsq - 2*radsq*(st1*st2*c12 + ct1*ct2)
#	print "hi: ",(st1*st2*c12 - ct1*ct2)
#	print "drsq: ", drsq
	dr = np.sqrt(drsq)
	drsqrt = np.sqrt(dr)
#	print dr
	drinv = 1/dr
	fx += pair[0][2]*pair[1][2]*drinv
	dvdr = pair[0][2]*pair[1][2]/drsq*(-1)
	drdt1 = 0.5*drinv*(-2*radsq*(st2*c12*ct1 - ct2*st1))
	drdt2 = 0.5*drinv*(-2*radsq*(st1*c12*ct2 - ct1*st2))
	drdp1 = 0.5*drinv*(+2*radsq*(st1*st2*s12))
	drdp2 = -drdp1
	i1 = pair[0][3]
	i2 = pair[1][3]
	dfx[i1] = dfx[i1] + drdt1*dvdr
	dfx[i2] = dfx[i1] + drdt2*dvdr
	dfx[i1+npts] = dfx[i1+npts] + drdp1*dvdr
	dfx[i2+npts] = dfx[i2+npts] + drdp2*dvdr

    return fx, dfx
    

    
    

# Main

def main():
    x0,q = genPts(npts)
    q[0] = 2
    q[1] = 2
    steps=int(sys.argv[1])
    myfile="hi.lammpstrj"
    g=open(myfile, "w+")
    g.close()
#    print x
#    x0 = np.array([0,np.pi*0.5,np.pi*0.25,np.pi*0.5])
#    xinit = x0
#    printRes(x0,q)
    writeLammps(x0, q, myfile)
    fx, dfx = coulomb(x0, q)
    print "fx: ", fx
    print "dfx: ", dfx
    for i in dfx:
	print i
    args=(q,)
    minimizer_kwargs = {"method": "BFGS", "jac":True, "args": args}
#    res = optimize.minimize(coulomb, x0, args, jac=True, options={'disp': True})
    res = optimize.basinhopping(coulomb, x0, minimizer_kwargs=minimizer_kwargs,niter=steps)
    print "global minimum: ", res.x, res.fun
#    x1 = res.x
#    fx, dfx = coulomb(x1, q)
#    print "fx: ", fx
#    print "dfx: ", dfx
#    printRes(res.x,q)
    writeLammps(res.x, q, myfile)
#    

if __name__ == "__main__":
    main()


