"""
Tests that the code reproduces the published curves of Martinez et al. (2016),
Fig. 2(a), at their parameters.

The known answer here is a published figure, read by pixel analysis to about
0.002 in nu (see schwinger/replication.py). The tolerance is that reading
uncertainty, not something chosen to make the tests pass.
"""
import numpy as np
import pytest

pytest.importorskip("qiskit")

from schwinger import replication as R

TOL = R.READ_UNCERTAINTY


def test_their_pulse_order_gives_the_same_density_as_ours():
    """
    Their step applies the hopping first and the diagonal terms last; ours the
    reverse. From the bare vacuum the two orders give the same nu at every
    step, built here two independent ways (Qiskit circuit, scipy expm).
    """
    assert np.allclose(R.martinez_pulse_order_density(), R.trotter_density(), atol=1e-12)


def test_reproduces_their_ideal_curve():
    """Exact nu(t) against the 33 points read off their green curve."""
    ideal = np.array(R.FIG2A_IDEAL)
    assert np.max(np.abs(R.exact_density(ideal[:, 0]) - ideal[:, 1])) < TOL


def test_reproduces_their_discretization_curve():
    """First-order Trotter steps of pi/8 against their red triangles."""
    assert np.max(np.abs(R.trotter_density() - np.array(R.FIG2A_TROTTER))) < TOL


def test_reproduces_their_error_model_curve():
    """Their fitted dephasing, p = 0.038 per qubit per step, against their squares."""
    assert np.max(np.abs(R.error_model_density() - np.array(R.FIG2A_MODEL))) < TOL


def test_a_different_step_length_does_not_reproduce_it():
    """
    The match identifies the step: pi/8 fits, while a step 10 % longer does not.
    Without this the replication tests could pass for any nearby step.
    """
    wrong = R.trotter_density(step=1.1 * R.STEP)
    assert np.max(np.abs(wrong - np.array(R.FIG2A_TROTTER))) > 10 * TOL


def test_their_dephasing_model_is_what_moves_the_last_point():
    """At p = 0 the model reduces to the Trotter curve, which misses their squares."""
    assert np.allclose(R.error_model_density(p=0.0), R.trotter_density(), atol=1e-12)
    assert abs(R.trotter_density()[-1] - R.FIG2A_MODEL[-1]) > 20 * TOL
