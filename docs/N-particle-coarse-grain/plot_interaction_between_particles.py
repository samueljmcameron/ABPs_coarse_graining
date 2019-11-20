import numpy as np
import matplotlib.pyplot as plt

def b_0(x):

    return np.pi/20*(36/2**(5/3)+8/x**10-20/x**4-10*x**2)

def b_2(x):

    return -np.pi/8*(18/2**(4/3)+2/x**8-8/x**2-x**4)


fig,ax = plt.subplots()
fig.set_size_inches(3.37,3.37)


xs = np.linspace(0.6,3,num=10000,endpoint=True)

ax.plot(xs,b_0(xs),label=r"$b_0/\sigma^2$")
ax.plot(xs,b_2(xs),label=r"$b_2/\sigma^4$")

ax.set_ylim(-10,10)
ax.set_xlim(0,3.5)

ax.legend(frameon=False,loc="upper right")


ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')

ax.set_xlabel(r"$a/\sigma$",x=1,labelpad=-1)
ax.set_xticks([1,2,3])
ax.set_yticks([-10,0,10])

fig.savefig("Figures/WCA_naive_b_0_b_2.pdf")
