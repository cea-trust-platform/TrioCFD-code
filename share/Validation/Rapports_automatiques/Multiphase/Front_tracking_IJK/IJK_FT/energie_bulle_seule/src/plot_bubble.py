import numpy as np
import matplotlib.pyplot as plt

DIST = 0.00008

def get_normal_vect(x, y):
    vect = np.array([x[1] - x[0], y[1] - y[0]])
    norm = np.sum(vect**2.)**0.5
    vect_norm = np.array([vect[1], -vect[0]]) / norm
    return vect_norm

if __name__ == "__main__":
    coo = np.loadtxt('bulle_err.txt')
    nodes = np.loadtxt('bulle_err_c.txt', dtype=int)

    plt.figure()
    for i in range(nodes.shape[0]):
        vect_norm = get_normal_vect(coo[nodes[i], 0], coo[nodes[i], 1])
        alea = np.random.uniform(0.5, 1., 1)
        vect_norm *= DIST*alea
        c, = plt.plot(coo[nodes[i], 0] + vect_norm[0], coo[nodes[i], 1] + vect_norm[1])
        plt.plot(coo[nodes[i], 0], coo[nodes[i], 1], c=c.get_color())

    plt.gca().set_aspect('equal')
    plt.show()
