# ZK Range Proof Implementation

A Python implementation of zero-knowledge range proofs, allowing a prover to convince a verifier that a secret number lies within a specific range without revealing the number itself.

- `run.py`: Main script demonstrating the range proof protocol
- `prover.py` and `verifier.py`: Prover and Verifier class implementations
- `utils.py`: Helper functions for elliptic curve operations, vector arithmetic and such.

This repo is based on learnings from the free-to-read [Rareskills ZK Book (Bulletproofs)](https://www.rareskills.io/zk-book).

## Prerequisites

- Python 3.8 or higher
- (Recommended) `uv` Python package manager

## Installation

1. Clone the repository:

```bash
git clone https://github.com/nkrishang/zk-range-proof.git
cd zk-range-proof
```

2. Create and activate a virtual environment:

```bash
# Using uv (recommended)
uv venv
```

3. Install dependencies

```bash
# Using uv (recommended)
uv pip install libnum py-ecc
```

## Usage

```bash
# Using uv (recommended)
uv run run.py
```

This will create a zero knowledge proof that the provided `v` and `n` are such that `v < 2^n` and verify the proof. It prints "accepted" if the verification is successful.
