# 🔐 Dual-Noise Encryption for DNA Storage

This repository implements a **modulation-based encryption scheme** that combines:

- Artificially injected errors (controllable, security-enhancing)
- Natural channel noise (random, modeled as sequencing errors)

It provides a complete pipeline including **key generation**, **error-injected encoding**, **decryption** (both legitimate and under attack), and **evaluation**.

---

## 📁 Project Structure

```
dual_noise_encryption/
│
├── encoding/                         # Encode plaintext into noisy DNA strands
│   └── encode-to-base-errorfocus.py   # Inject artificial + natural channel noise
│
├── decoding/                         # Decode noisy DNA reads
│   ├── clusteranalysis.py              # Clustering and read preprocessing
│   ├── withkey_decryption.py          # Legitimate decryption using correct key
│   ├── MSA_decryption_direct.py       # Chosen-plaintext attack (Direct MSA)
│   ├── MSA_decryption_direct_double.py
│   ├── MSAinfer-decryption.py         # Inferred-key attack via MSA
│   ├── MSAinfer-decryption_double.py
│   ├── cluster_double.py              # Attack-aware clustering + decryption
│   └── mafft-win/                     # MAFFT executable for sequence alignment
│
├── key_generation/                   # Generate constrained modulation keys
│   ├── generation_main.py              # Biochemically valid + misleading key generation
│   └── generated_key_2_5.txt           # Example key file
```

---

## 🔧 Functional Overview

### 1. Key Generation (`key_generation/`)

- Generates **modulation keys**: binary sequences that satisfy biochemical constraints (e.g., max homopolymer length).
- Also creates **misleading keys** (decoys) with fixed Hamming distance to the real key.

**Script**: `generation_main.py`

---

### 2. Encoding with Dual Noise (`encoding/`)

- Converts plaintext into binary.
- Encodes binary with modulation keys into DNA bases.
- **Adds two types of noise**:
  - **Artificially injected noise** (structured, attacker-controlled)
  - **Channel noise** (random substitution, insertion, deletion)

**Script**: `encode-to-base-errorfocus.py`

---

### 3. Decoding & Attack Simulation (`decoding/`)

- **Legitimate decryption** (`withkey_decryption.py`):
  - Assumes the correct modulation key is available.
  - Performs clustering, alignment, and demodulation.

- **Adversarial decryption**:
  - **Direct MSA attack** (`MSA_decryption_direct.py`, `_double.py`):
    - No clustering; assumes access to chosen plaintext.
  - **Inferred-key attack** (`MSAinfer-decryption.py`, `_double.py`):
    - Infers the modulation key from the noisy read pool.

- **Clustering tool** (`clusteranalysis.py`):  
  - Groups noisy reads for further processing.

- **MAFFT Integration**:
  - All MSA-based decryption uses the `mafft-win/mafft.bat` executable.

---

## 🧪 Requirements

- Python 3.7+
- [Biopython](https://biopython.org/) (`pip install biopython`)
- [Crypto](https://pypi.org/project/pycryptodome/) for secure random generation
- `matplotlib`, `numpy`, `pandas`
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (Windows batch file included)

---

## 📝 Usage Flow

1. **Generate keys**  
   `python generation_main.py`  
   → Outputs `generated_key_2_5.txt`

2. **Encode message with noise**  
   `python encode-to-base-errorfocus.py --error_rate 0.05 --errorrate 0.1`  
   → Simulates strand synthesis/sequencing

3. **Decrypt (or attack)**  
   Use:
   - `withkey_decryption.py` for legit decoding
   - `MSA_decryption_direct.py` / `MSAinfer-decryption.py` for adversarial

4. **Evaluate**  
   Output includes:
   - Bit error rate (BER)
   - Character error rate (CER)
   - Key inference accuracy

---

## 📌 Notes

- File paths are currently Windows-style (`D:\...`) — adjust if on Linux/macOS.
- The `_double.py` files simulate **dual clustering** attack scenarios.
- `mafft.bat` must be accessible in `decoding/mafft-win`.

---

## 🧠 Citation & Acknowledgement

If you use this code in your work, please cite the related publication (TBD).  
MAFFT is developed by Katoh et al., and should be cited accordingly.
