[Tour](README.md) · *1 of 9* · [Next: from a line of space to four qubits →](02-four-sites-four-qubits.md)

# 1. Why simulate empty space?

Put enough energy into one small region and an electron and its antiparticle,
the positron, appear. This is Einstein's $E = mc^2$ read backwards: energy
turning into mass. Particle detectors record it routinely, and the theory
behind it, quantum electrodynamics (QED), predicts it in detail.

Following that process as it happens, moment by moment, is another matter.
This project does it for the simplest version of QED there is, on a simulated
quantum computer, and then asks how much of the computer's answer can be
trusted.

## One dimension instead of three

The Schwinger model is QED with one dimension of space instead of three
([Schwinger 1962](https://doi.org/10.1103/PhysRev.128.2425)). It keeps the two
features this project cares about.

*Pairs.* Particles and antiparticles can be created and destroyed.

*Confinement.* In three dimensions the field lines of a charge spread out
over a sphere, so the force falls off as $1/r^2$. On a line they have nowhere
to spread. The field between a plus and a minus charge stays the same however
far apart they are, so the energy stored in it grows in proportion to the
distance ([Coleman, Jackiw and Susskind 1975](https://doi.org/10.1016/0003-4916(75)90212-2)).
Pull the two charges apart and the cost keeps rising until it is cheaper to
make a new pair out of the field. The new charges screen the old ones and the
"string" of field snaps.

Quarks behave the same way in three dimensions: no experiment has ever pulled
a single quark out of a proton. That is why this model is the standard test
bed for the theory of quarks, quantum chromodynamics (QCD). The particles here
are electrons and positrons, not quarks. What they share with quarks is the
confinement.

## Why an ordinary computer struggles

There are two separate problems.

*Memory.* Divide the line into $N$ sites, each of which either holds a
particle or does not. The quantum state is a list of one complex number (an
amplitude) for every possible arrangement, because the system can be in all
of them at once. That is $2^N$ numbers. Four sites need 16. Fifty sites need
about $10^{15}$. Every site you add doubles the list.

*Sampling fails in real time.* The method that computes proton masses from
QCD, lattice Monte Carlo, samples field configurations the way a pollster
samples a population. It works because it replaces real time by imaginary
time ($t \to -i\tau$), which turns the weight of every configuration into a
positive probability. Ask about real time instead and the weights become
complex phases that mostly cancel, and the statistical error grows
exponentially with time and system size. This is the sign problem
([Troyer and Wiese 2005](https://doi.org/10.1103/PhysRevLett.94.170201)), and it
is why watching pairs form in real time is out of reach of the methods that
work so well for static properties.

## Why a quantum computer, and why that is not the end of it

A quantum computer holds the state in its qubits directly: $N$ qubits carry a
state of $2^N$ amplitudes. It then applies the time evolution as a sequence
of gates. Feynman proposed simulating quantum physics this way in 1982
([Feynman 1982](https://doi.org/10.1007/BF02650179)), and Lloyd showed in 1996
that it works for any Hamiltonian built from local pieces, with an error you
can make as small as you like by cutting time into short enough steps
([Lloyd 1996](https://doi.org/10.1126/science.273.5278.1073)).

Today's machines make mistakes, though, and every answer one gives carries
three kinds of error.

1. *Time steps.* The evolution is applied in small chunks, which is exact
   only in the limit of infinitely many chunks. This is Trotter error, and it
   would be there even on a perfect machine.
2. *Shots.* Measuring a qubit gives one random bit. An average needs many
   runs of the circuit (shots), and its uncertainty falls only as one over the
   square root of their number.
3. *The hardware.* Gates are imperfect and readout sometimes misreports a
   bit.

The three pull against each other. Smaller time steps cut the first error but
need more gates, which raises the third.

## The experiment this project rebuilds

In 2016 a group in Innsbruck ran this model on four trapped-ion qubits and
watched electron-positron pairs appear out of the vacuum
([Martinez et al. 2016, Nature 534, 516](https://doi.org/10.1038/nature18318)).
It was the first digital quantum simulation of a lattice gauge theory. Their
measured points drift away from the exact curve as time goes on. How much of
that gap does each kind of error cause?

On real hardware that cannot be untangled: the errors arrive together and the
exact answer is unknown. So this project rebuilds the experiment on a
simulated quantum computer, where the exact answer is known at every stage,
switches the errors on one at a time, and measures what each one does. It is
after two numbers and one test.

- *The crossover.* At what time step does the error from chopping time fall
  below the shot noise? Past that point, smaller steps buy nothing you can see.
- *The failure threshold.* How bad can the gates get before the
  pair-production signal is lost in the noise?
- *The warning light.* Can the machine tell you it has failed, using only its
  own data, without knowing the right answer?

## Why only four sites

Four sites give a 16-by-16 Hamiltonian matrix, small enough that the exact
answer is available at every stage. It is also the size of the 2016
experiment, so the results line up with real data. And the noisy simulations
are cheap enough to repeat across dozens of error rates and step sizes.

At four sites a laptop beats any quantum computer, so nothing here is a claim
of quantum advantage. The point is the method: these are the checks you would
want once the problem is too big to check. [Chapter 3](03-why-pairs-appear.md)
shows that four sites already get the first swing of pair production right,
by solving the same model exactly on up to twelve sites.

---

*In the thesis:* Section 1, Introduction ([PDF](../../thesis/main.pdf)).
*Deeper:* [`docs/derivation.md`](../derivation.md), sections 0 and 1.

[Tour](README.md) · [Next: from a line of space to four qubits →](02-four-sites-four-qubits.md)
