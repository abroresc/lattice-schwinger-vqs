# Setup, and getting this onto GitHub

## 1. Get it running on your laptop, once

You need Python 3.11 or newer, VS Code, and git.

1. **Python.** python.org, or `brew install python` on a Mac. Check with
   `python3 --version`.
2. **VS Code.** code.visualstudio.com. Then install two extensions from the
   Extensions panel: **Python** and **Jupyter**, both by Microsoft.
3. **Open the project.** File > Open Folder, pick `schwinger-n4`.
4. **Make a virtual environment.** In the VS Code terminal
   (Terminal > New Terminal):

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate        # Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   ```

   Then bottom-right of VS Code, or Cmd/Ctrl+Shift+P > "Python: Select
   Interpreter", pick the one inside `.venv`.

5. **Check it works.**

   ```bash
   python -m pytest tests/ -q
   ```

   You should get 394 passed. If you do, everything below will work.

## 2. Day to day

```bash
python -m pytest tests/ -q                    # after any change to src/
python scripts/03_schwinger_spectrum.py       # the Gate 2 artifact
python scripts/05_qiskit_trotter.py           # circuits, Trotter, shot noise
```

Open `notebooks/schwinger_walkthrough.ipynb` in VS Code and it runs inline. That
notebook is the thing to show people.

**Why VS Code and not Jupyter Lab on its own.** This is a package with several
modules and a test suite, not a single scratch file. VS Code edits the package,
runs the tests in a terminal, opens the notebook natively, and does git, all in
one window. Jupyter Lab does the notebook well and everything else badly. You
lose nothing by using VS Code: `.ipynb` files open there with the same cell
interface.

**Keep the exploration in notebooks and the real code in `src/`.** Anything you
want to keep and test goes into a module; the notebook imports it. Code that
lives only in notebook cells is how projects become undefendable.

## 3. Putting it on GitHub

Git history is already initialised with a first commit, so this is three steps.

1. Make an empty repo at github.com/new. Name it `schwinger-n4`. **Do not** tick
   "Add a README", "Add .gitignore", or "Choose a license", since this repo
   already has all three and it will refuse to push.
2. Decide public or private. Nothing here is sensitive.
   You can start private and flip it public later from Settings.
3. In the VS Code terminal, from inside `schwinger-n4`:

   ```bash
   git remote add origin https://github.com/<your-username>/schwinger-n4.git
   git branch -M main
   git push -u origin main
   ```

   It will ask you to sign in. Use the browser flow or a personal access token.

After that, `git add -A && git commit -m "..." && git push` whenever you finish
something.

**Set the author name first** if you have not used git before, so the commits
are yours:

```bash
git config --global user.name "Abror Niyazmetov"
git config --global user.email "you@example.com"
```

Then, for the existing commit only, rewrite the author:

```bash
git commit --amend --reset-author --no-edit
```

