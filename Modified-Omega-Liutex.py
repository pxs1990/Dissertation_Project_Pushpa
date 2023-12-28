import numpy as np

# Function for rotation matrix
def rotation(u, v):
    a = np.cross(u, v)
    aa = np.linalg.norm(a)
    t = np.dot(u, v)
    alpha = np.arccos(t / (np.linalg.norm(u) * np.linalg.norm(v)))

    if aa != 0:
        c = aa
        s = np.sin(alpha) * aa
        r = np.array([
            [c + a[0] * a[0] * (1 - c),      a[0] * a[1] * (1 - c) - a[2] * s, a[0] * a[2] * (1 - c) + a[1] * s],
            [a[1] * a[0] * (1 - c) + a[2] * s, c + a[1] * a[1] * (1 - c),      a[1] * a[2] * (1 - c) - a[0] * s],
            [a[2] * a[0] * (1 - c) - a[1] * s, a[2] * a[1] * (1 - c) + a[0] * s, c + a[2] * a[2] * (1 - c)]
        ])
    else:
        r = np.identity(3)

    return r

# Example usage:
u = np.array([1.0, 0.0, 0.0])
v = np.array([0.0, 1.0, 0.0])

rotation_matrix = rotation(u, v)
print(rotation_matrix)
