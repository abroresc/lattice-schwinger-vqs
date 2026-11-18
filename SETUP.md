# Setup

Everything here runs on a laptop. Four sites make a 16-dimensional problem, so
nothing needs a cluster or a quantum computer; the slowest script takes about
five minutes.

## 1. Get it running

You need Python 3.11 or newer and git.

```bash
git clone https://github.com/abroresc/lattice-schwinger-vqs.git
cd lattice-schwinger-vqs
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

`requirements.txt` gives minimum versions. The exact versions every number in
the thesis came from are pinned in `requirements-lock.txt`; install from that
file instead if a residual ever differs from the published one.

Check that it works:

```bash
python -m pytest tests/ -q
```

You should see `496 passed` in about fifteen seconds. If you do, everything
below will work.

> **If you use Anaconda,** a bare `python` may be Anaconda's, with different
> library versions, even inside the virtual environment. Call
> `.venv/bin/python` explicitly to be sure.

## 2. Exploring

```bash
python scripts/03_schwinger_spectrum.py      # the Hamiltonian and its checks
python scripts/04_pair_production_preview.py # pairs out of the vacuum
python scripts/05_qiskit_trotter.py          # circuits, Trotter error, shots
python scripts/06_noise_threshold.py         # gate error and the failure threshold (slow)
bash scripts/verify_all.sh                   # everything, rewriting VERIFICATION.txt
```

The README has a table of every script and what it shows. The notebook,
`notebooks/schwinger_walkthrough.ipynb`, opens in Jupyter or in VS Code with
the Python and Jupyter extensions.

## 3. The thesis

```bash
cd thesis && make        # builds main.pdf and intro_methods.pdf
```

This needs a LaTeX installation with `latexmk` (BasicTeX on a Mac is enough).
Without one, upload `thesis/` and `figures/` to Overleaf and press Recompile.
The thesis includes the figure PDFs straight from `figures/`, so rerunning a
script updates the thesis on the next build.

## 4. If you change something

- Run `python -m pytest tests/ -q` after any change to `src/`.
- New code gets a test against an answer known independently: a hand-written
  matrix, an algebraic identity, a closed-form result, or a second
  construction that shares no code with the first.
- Exploration can live in a notebook; anything worth keeping moves into
  `src/` with a test.
- Read [`CONVENTIONS.md`](CONVENTIONS.md) before touching site order, signs or
  Qiskit qubit indices. Those are the conventions that give plausible wrong
  answers without any error message.
