#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
from subprocess import Popen, PIPE
import shlex

def get_command_output(command):
    p = Popen(shlex.split(command), stdin = PIPE, stdout = PIPE,
        stderr = PIPE)
    (stdout, stderr) = p.communicate()
    strout = bytes.decode(stdout)
    if(p.returncode == 0):
        str1D = strout.split()
        r      = np.array(list(map(float, str1D[0::4])))
        numRes = np.array(list(map(float, str1D[1::4])))
        anaRes = np.array(list(map(float, str1D[2::4])))
        return (r, numRes, anaRes)
    else:
        strerr = bytes.decode(stderr)
        print('ERROR: the coammnd <' + command + '> goes wrong!')
        print(strout + strerr)
        exit()

if len(sys.argv) > 1:
    ie = float(sys.argv[1])
else:
    ie = 4

(rGup, numGup, anaGup) = get_command_output('./Guptasarma')
errGup = np.abs((numGup - anaGup)/anaGup)

(rOga, numOga, anaOga) = get_command_output('./Ogata')
errOga = np.abs((numOga - anaOga)/anaOga)

if(ie == 4):
    equ = r'$ \int_0^{\infty} e^{ - c \lambda} J_0 (r \lambda) d \lambda ' \
        + r'= \frac{1}{ \sqrt{ c^2 + r^2 } } $'
elif(ie == 5):
    equ = r'$ \int_0^{\infty} \lambda e^{ - c \lambda^2} J_0 (r \lambda) ' \
        + r'd \lambda = \frac{1}{2c} e^{ - \frac{r^2}{4 c} } $'
elif(ie == 6):
    equ = r'$ \int_0^{\infty} \lambda e^{ - c \lambda} J_0 (r \lambda) ' \
        + r'd \lambda = \frac{c}{ \left( c^2 + r^2 \right)^{3/2} } $'
elif(ie == 7):
    equ = r'$ \int_0^{\infty} \left[ \lambda e^{ - c \lambda} + \alpha ' \
        + r'\lambda^2 e^{ - c \lambda^2} \right] J_1 (r \lambda) d \lambda ' \
        + r'= \frac{r}{ \left( c^2 + r^2 \right)^{3/2} } + \alpha ' \
        + r'\frac{ r e^{ - \left( r^2 / 4c \right) } }{4 c^2} $'
elif(ie == 8):
    equ = r'$ \int_0^{\infty} \lambda e^{ - c \lambda} J_1 (r \lambda) ' \
        + r'd \lambda = \frac{r}{ \left( c^2 + r^2 \right)^{3/2} } $'
elif(ie == 9):
    equ = r'$ \int_0^{\infty} \lambda^2 e^{ - c \lambda^2} J_1 (r \lambda) ' \
        + r'd \lambda = \frac{ r e^{ - \left( r^2 / 4c \right) } }{4 c^2} $'
elif(ie == 10):
    equ = r'$ \int_0^{\infty} e^{ - c \lambda} J_1 (r \lambda) d \lambda ' \
        + r'= \frac{ \sqrt{r^2 + c^2} - c }{ r \sqrt{r^2 + c^2} } $'
else:
    print('ERROR: no the example number <%d>'%(ie))
    exit()

plt.loglog(rGup, errGup, label = 'Guptasarma')
plt.loglog(rOga, errOga, label = 'Ogata')
plt.xlim((min(np.min(rGup), np.min(rOga)), max(np.max(rGup), np.max(rOga))))
plt.xlabel('r')
plt.ylabel('Relative error')
plt.legend(loc = 'upper left')

if 'equ' in vars():
    if 0:
        plt.title(equ, fontsize = 20)
    else:
        plt.text(7.0e-7, 5.0e-0, equ, {'fontsize': 20},
            verticalalignment = 'bottom', horizontalalignment = 'left')

if 1:
    plt.show()
else:
    plt.savefig('figures/example-%d.png'%(ie), dpi = 150)

