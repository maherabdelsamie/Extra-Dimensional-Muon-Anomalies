import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

cp.random.seed(123)

def expm_gpu(matrix):
    if cp.isscalar(matrix):
        return cp.exp(-1j * matrix)
    else:
        w, v = cp.linalg.eigh(matrix)
        exp_diag = cp.diag(cp.exp(-1j * w))
        return cp.matmul(cp.matmul(v, exp_diag), cp.linalg.inv(v))

def kinetic_energy(N):
    return cp.diag(cp.arange(N, dtype=cp.complex128))

def potential_energy(N):
    return cp.diag(cp.arange(N, dtype=cp.complex128)**2)

def gf_fluctuations(α, N):
    return cp.random.normal(scale=α, size=(N,N)).astype(cp.complex128)

def gt_transitions(β, N):
    return β * cp.diag(cp.arange(N+1, dtype=cp.complex128)[:N])

def gc_correlations(γ, psi2):
    return cp.real(γ * cp.outer(psi2, psi2))[:N,:N].astype(cp.complex128)

def run_simulation(α, β, γ, dt, Nt, N):
    mu_values = cp.zeros(Nt, dtype=cp.complex128)
    psi = cp.zeros(N, dtype=cp.complex128)
    psi[0] = 1
    psi2 = cp.random.randn(N) + 1j*cp.random.randn(N)
    psi2 /= cp.linalg.norm(psi2)
    for t in range(Nt):
        H = kinetic_energy(N) + potential_energy(N)
        H += gf_fluctuations(α, N)
        H += gt_transitions(β, N)
        H += gc_correlations(γ, psi2)
        psi = expm_gpu(-1j * H).dot(psi)

        # Normalize psi
        psi /= cp.linalg.norm(psi)

        mu_op = cp.diag(cp.arange(N, dtype=cp.complex128))
        mu_values[t] = cp.vdot(psi, mu_op.dot(psi)).real
    return mu_values, H

RANGE_MIN = 1.16589671e-8
RANGE_MAX = 1.16594439e-8

def main_simulation(alpha_values, beta_values, gamma_values, dt, Nt, N):
    anomalous_moments = []
    time_evolution_data = []
    hamiltonian_data = []
    total_iterations = len(alpha_values) * len(beta_values) * len(gamma_values)
    with tqdm(total=total_iterations) as pbar:
        for alpha in alpha_values:
            for beta in beta_values:
                for gamma in gamma_values:
                    mu_values, H = run_simulation(alpha, beta, gamma, dt, Nt, N)
                    time_evolution_data.append(mu_values)
                    hamiltonian_data.append(H)
                    mu_freq = cp.fft.fft(mu_values)
                    a_mu = mu_freq[1].real
                    anomalous_moments.append((alpha, beta, gamma, a_mu))
                    pbar.update(1)
                    if RANGE_MIN <= a_mu <= RANGE_MAX:
                        print(f"Computed Anomalous magnetic moment for α = {alpha}, β = {beta}, γ = {gamma}: {a_mu} (Within Range)")
                    else:
                        print(f"Computed Anomalous magnetic moment for α = {alpha}, β = {beta}, γ = {gamma}: {a_mu} (Not Within Range)")
    return anomalous_moments, time_evolution_data, hamiltonian_data

if __name__ == '__main__':
    dt = 1e-2
    Nt = int(1e4)
    N = 50

    alpha_values = cp.array([1.6237767391887243e-10])
    beta_values = cp.array([3.4e-9])
    gamma_values = cp.array([2.86057797e-9])
    anomalous_moments, time_evolution_data, hamiltonian_data = main_simulation(alpha_values, beta_values, gamma_values, dt, Nt, N)

# PLOTTING RESULTS
# 1. Time Evolution of Muon Observable
plt.figure(figsize=(10, 6))
for alpha, mu_values in zip(alpha_values.get(), time_evolution_data):
    plt.plot(mu_values.get(), label=f"α = {alpha:.2e}")
plt.title("Time Evolution of Muon Observable")
plt.xlabel("Time Steps")
plt.ylabel("Muon Observable Value")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("time_evolution.png")
plt.show()

# 2. Anomalous Magnetic Moment vs. α
plt.figure(figsize=(10, 6))
plt.plot(alpha_values.get(), [am[3].get() for am in anomalous_moments], 'o-')
plt.axhline(y=RANGE_MIN, color='r', linestyle='--', label="Experimental Min Value")
plt.axhline(y=RANGE_MAX, color='b', linestyle='--', label="Experimental Max Value")
plt.title("Anomalous Magnetic Moment vs. α")
plt.xlabel("α")
plt.ylabel("Anomalous Magnetic Moment")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("anomalous_vs_alpha.png")
plt.show()

# 3. Fourier Transform of Magnetic Moment
plt.figure(figsize=(10, 6))
for mu_values in time_evolution_data:
    frequencies = cp.fft.fftfreq(mu_values.size, d=dt)
    μ_freq = cp.fft.fft(mu_values, axis=0)
    plt.plot(frequencies.get(), cp.abs(μ_freq).get())
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.title('Fourier Transform of Magnetic Moment')
plt.tight_layout()
plt.grid(True)
plt.show()

# 4. Probability Density Evolution
plt.figure(figsize=(10, 6))
for alpha, mu_values in zip(alpha_values.get(), time_evolution_data):
    plt.plot((mu_values**2).get(), label=f"α = {alpha:.2e}")
plt.title("Probability Density Evolution")
plt.xlabel("Time Steps")
plt.ylabel("Probability Density")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("probability_density_evolution.png")
plt.show()

# 5. Hamiltonian Components Visualization
plt.figure(figsize=(10, 6))
alpha_sample = alpha_values[0].get()
beta_sample = beta_values[0].get().item()  # Extract scalar value
gamma_sample = gamma_values[0].get().item()  # Extract scalar value
plt.plot(kinetic_energy(N).diagonal().get(), label="Kinetic Energy")
plt.plot(potential_energy(N).diagonal().get(), label="Potential Energy")
plt.plot(gf_fluctuations(alpha_sample, N).diagonal().get(), label="GF Fluctuations")
plt.plot(gt_transitions(beta_sample, N).diagonal().get(), label="GT Transitions")
plt.plot(gc_correlations(gamma_sample, cp.zeros(N, dtype=cp.complex128)).diagonal().get(), label="GC Correlations")
plt.title("Hamiltonian Components for α = {}".format(alpha_sample))
plt.xlabel("State Index")
plt.ylabel("Component Value")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("hamiltonian_components.png")
plt.show()


# Fourier Transform of Contributions
T = kinetic_energy(N)
V = potential_energy(N)
gf = gf_fluctuations(alpha_sample, N)
gt = gt_transitions(beta_sample, N)
gc = gc_correlations(gamma_sample, cp.zeros(N, dtype=cp.complex128))
    
# Compute the Fourier Transform of each component
T_freq = cp.fft.fft(T.diagonal())
V_freq = cp.fft.fft(V.diagonal())
gf_freq = cp.fft.fft(gf.diagonal())
gt_freq = cp.fft.fft(gt.diagonal())
gc_freq = cp.fft.fft(gc.diagonal())
    
# Plot the Fourier Transform of each component
plt.figure(figsize=(10, 6))
plt.plot(cp.abs(T_freq).get(), label="Kinetic Energy")
plt.plot(cp.abs(V_freq).get(), label="Potential Energy")
plt.plot(cp.abs(gf_freq).get(), label="GF Fluctuations")
plt.plot(cp.abs(gt_freq).get(), label="GT Transitions")
plt.plot(cp.abs(gc_freq).get(), label="GC Correlations")
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.legend()
plt.title("Fourier Transform of Contributions")
plt.show()
