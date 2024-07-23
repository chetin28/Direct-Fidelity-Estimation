import numpy as np
from qiskit import *
from qutip import *
from qiskit.opflow import I, X, Y, Z, PauliSumOp
import qiskit.quantum_info as qi
import math
from qiskit.tools import job_monitor

IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q-education', group='mid-east-tech-un-1', project='2300343-Intro-Computational-Methods')
backend = provider.get_backend("ibmq_belem")


n = 5
d = 2**n

q = QuantumRegister(5)
c = ClassicalRegister(5)
qc = QuantumCircuit(q,c)

#pure state rho, its not basis(d,0) so that qutip can know we are working on 5 qubit system
rho=tensor(basis(2,0),basis(2,0),basis(2,0),basis(2,0),basis(2,0))

epsilon = 0.5
delta = 0.5

l = round(1/(epsilon**2 * delta))

def matrixinit(X,xi): # X and xi are lists
    if len(xi) == 5:
        for i in range(len(xi)):
            k = xi[i]
            kk = k[1]
            X.append(exp_vals[kk])
            
        return X
    else:
        for i in range(len(xi)):
            k = xi[i]
            kk = k[1]
            kkk = k[0]
            ii = 0
            if i == 0:
                ii = i
            else:
                ii = len(X)
            while ii<kkk+1:
                if ii != kkk:
                    X.append(0)
                    ii = ii + 1
                if ii == kkk:
                    X.append(exp_vals[kk])
                    break
        if len(X)!=5:
            while len(X)!=5:
                X.append(0)  
        return X

omega =[]
omega_tilde = []
Aij = []

omega_tilde_i = 0
pauliset = [qeye(2),sigmax(),sigmay(),sigmaz()]

while len(omega_tilde) < l:
    chi_rho_k = 0
    while chi_rho_k == 0:
        rand=np.random.randint(0,4,size=n)
        #construct W_k
        Wk = tensor(pauliset[rand[0]],pauliset[rand[1]],pauliset[rand[2]],pauliset[rand[3]],pauliset[rand[4]])
        chi_rho = Wk*rho/math.sqrt(d)
        chi_rho_k = chi_rho.tr()
           #implementation on the circuit for sigma
            # m: number of required copies of sigma
    m = math.ceil((2*math.log(2/delta)) / (chi_rho_k**2 * l * epsilon**2 * d))
    
    for sth in range(m): # define m numbers of circuits
        chi_rho_k = 0
        while chi_rho_k == 0:
            rand=np.random.randint(0,4,size=5)
            Wk = tensor(pauliset[rand[0]],pauliset[rand[1]],pauliset[rand[2]],pauliset[rand[3]],pauliset[rand[4]])
            chi_rho = Wk*rho/math.sqrt(d)
            chi_rho_k = chi_rho.tr()
        W_k = []
        for i in range(5):
            if rand[i] == 0: W_k += "I"
            if rand[i] == 1: W_k += "X"
            if rand[i] == 2: W_k += "Y"
            if rand[i] == 3: W_k += "Z"
        print("chosen W_k="+str(W_k))
        qc.initialize("00000",qc.qubits)
        for i in range(5):
            if rand[i] == 0: qc.id(q[i])
            if rand[i] == 1: qc.x(q[i])
            if rand[i] == 2: qc.y(q[i])
            if rand[i] == 3: qc.z(q[i])
        qc.measure(q,c)
        
        #create job
        trans = transpile(qc , backend, optimization_level=2)
        job = backend.run(trans)
        job_monitor(job,interval=2)
        results = job.result()
        #simulation
        #job = execute(qc,Aer.get_backend("qasm_simulator"))
        #results=job.result().get_counts(qc)
        chi_sigma_real = results.get_counts()
        print("real:", chi_sigma_real)
        
        #reading and processing measurements
        values = list(chi_sigma_real.values())
        keys = list(chi_sigma_real.keys())
        total = 0
        exp_vals = []
        for i in range(len(values)):
            total = total + values[i]
        for i in range(len(values)):
            exp_vals.append(values[i]/total)
        label=[]
        for i in range(len(keys)):
            label.append(str(keys[i]))
        base_ten = []
        for i in range(len(label)):
            a = int(keys[i],2)
            base_ten.append(a)
        A,B,C,D,E = [],[],[],[],[]
        ai,bi,ci,di,ei = [],[],[],[],[]
        for i in range(len(exp_vals)):
            if base_ten[i]/20 > 1:
                eii = base_ten[i]%20
                ei.append([eii,i])
            if 1.3 > base_ten[i]/15 > 1:
                dii = base_ten[i]%15
                di.append([dii,i])
            if 1.45 > base_ten[i]/10 > 1:
                cii = base_ten[i]%10
                ci.append([cii,i])
            if 1.9 > base_ten[i]/5 > 1:
                bii = base_ten[i]%5
                bi.append([bii,i])
            if 5 > base_ten[i]:
                aii = base_ten[i]
                ai.append([aii,i])
        ai.sort(), bi.sort(), ci.sort(), di.sort(), ei.sort()
        
        prnta = matrixinit(A,ai)
        prntb = matrixinit(B,bi)
        prntc = matrixinit(C,ci)
        prntd = matrixinit(D,di)
        prnte = matrixinit(E,ei)
        mAij = [prnta,prntb,prntc,prntd,prnte] #matrix form of measurement results
        qobj = Qobj(mAij) #turning into quantum object
        trqobj = qobj.tr()
        Aij.append(trqobj)
        omega.append(trqobj/chi_rho_k)
        omega_tilde_i =  omega_tilde_i + trqobj / ( m * math.sqrt(d) * chi_rho_k  )
        print(trqobj)
    omega_tilde.append(omega_tilde_i)
print("A_ij=",Aij)


#infinite-precision estimator
xhi = 0
stop = 0
for i in range(len(omega)):
    xhi = xhi + omega[i]/(l*m*math.sqrt(d))

#print(omega_tilde[6])
xhi_tilde = 0
for i in range(l):
    xhi_tilde = xhi_tilde + omega_tilde[i]/(l*math.sqrt(d))
print(xhi_tilde)


Fmax = xhi_tilde + 2*epsilon
if Fmax > 1:
    Fmax = 1
Fmin = xhi_tilde - 2*epsilon
print("fidelity lies between [",Fmin,Fmax,"]")


#failure prob
if abs(xhi_tilde-xhi) >= epsilon:
    print("error occured")
    print("xhi: ", xhi)
    print("xhi_tilde: ",xhi_tilde)
else:
    print("completed with following difference and fixed error respectively")
    print(abs(xhi_tilde-xhi))
    print(epsilon)
