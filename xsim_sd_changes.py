import numpy as np

def simulate_sd_blocks(n, sd, k, seed=None):
    """Simulate k columns with blockwise standard deviations."""
    n = np.asarray(n, dtype=int)
    sd = np.asarray(sd, dtype=float)

    if n.ndim != 1:
        raise ValueError("n must be a 1d array")
    if sd.ndim != 1:
        raise ValueError("sd must be a 1d array")
    if n.size != sd.size:
        raise ValueError("n and sd must have the same size")
    if np.any(n < 0):
        raise ValueError("all elements of n must be nonnegative")
    if np.any(sd < 0.0):
        raise ValueError("all elements of sd must be nonnegative")

    rng = np.random.default_rng(seed)
    blocks = [rng.normal(loc=0.0, scale=sd[i], size=(n[i], k)) for i in range(n.size)]
    x = np.vstack(blocks)
    return x

def write_matrix_with_header(filename, x, n, sd, k, seed):
    """Write matrix to text file with parameter lines starting with #."""
    n = np.asarray(n, dtype=int)
    sd = np.asarray(sd, dtype=float)
    nrow, ncol = x.shape

    with open(filename, "w", encoding="ascii") as f:
        f.write(f"# {nrow} {ncol} nrow ncol\n")
        f.write("# distribution: normal\n")
        f.write(f"# n = {n.tolist()}\n")
        f.write(f"# sd = {sd.tolist()}\n")
        f.write("# mean = 0.0\n")
        f.write(f"# seed = {seed}\n")
        np.savetxt(f, x, fmt="%.8f")

def main():
    n = np.array([100, 120, 80], dtype=int)
    sd = np.array([1.0, 2.0, 0.5], dtype=float)
    k = 100
    seed = 12345
    filename = "sim_matrix.txt"

    x = simulate_sd_blocks(n, sd, k, seed=seed)
    write_matrix_with_header(filename, x, n, sd, k, seed)
    print(f"wrote {x.shape[0]} x {x.shape[1]} matrix to {filename}")

if __name__ == "__main__":
    main()
