#
#  Java Information Dynamics Toolkit (JIDT)
#  Copyright (C) 2026, Alexander Richter
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Checks that the module supplies a working toolkit.

The measures themselves are tested in java/unittests; what is checked here
is the wiring: that a JVM starts on the runtime supplied with the module,
that the classes in the jar are reachable by name, and that observations may
be given as Python data and results read back as such.
"""

import zipfile

import jpype
import numpy as np
import pytest

import jidt


@pytest.fixture
def data():
    """A random number generator, seeded so that tests do not vary."""
    return np.random.default_rng(0)


def test_jar_and_runtime_are_bundled():
    assert jidt.JAR.exists()
    assert (jidt.java_home() / "bin" / "java").exists()


def test_a_32_bit_python_is_rejected_with_a_clear_message(monkeypatch):
    monkeypatch.setattr(jidt.struct, "calcsize", lambda fmt: 4)
    with pytest.raises(RuntimeError, match="64-bit Python"):
        jidt.check_platform()
    with pytest.raises(RuntimeError, match="64-bit Python"):
        jidt.start()


def test_a_64_bit_python_passes_the_check():
    jidt.check_platform()


def test_classes_are_exposed_under_their_java_names():
    assert "MutualInformationCalculatorDiscrete" in dir(jidt)
    assert "TransferEntropyCalculatorKraskov" in dir(jidt)
    with pytest.raises(AttributeError):
        # The lookup itself is the assertion here:
        jidt.NoSuchCalculator  # noqa: B018


def test_lookup_by_full_java_name():
    assert jidt.jclass("infodynamics.utils.MatrixUtils") is not None


def test_observations_may_be_given_as_numpy_or_as_lists(data):
    """The same numbers must reach Java whichever container holds them."""
    x = data.integers(0, 2, 2000)
    y = data.integers(0, 2, 2000)

    calc = jidt.MutualInformationCalculatorDiscrete(2)
    calc.initialise()
    calc.addObservations(x, y)
    from_numpy = calc.computeAverageLocalOfObservations()

    calc.initialise()
    calc.addObservations(x.tolist(), y.tolist())

    assert calc.computeAverageLocalOfObservations() == from_numpy


def test_univariate_observations_match_a_hand_built_java_array(data):
    """Our conversion must agree with the Java arrays it stands in for."""
    source = data.normal(size=1000)
    destination = source + data.normal(size=1000) * 0.5

    ours = jidt.MutualInfoCalculatorMultiVariateKraskov1()
    ours.initialise(1, 1)
    ours.setObservations(source, destination)

    theirs = jidt.MutualInfoCalculatorMultiVariateKraskov1()
    theirs.initialise(1, 1)
    theirs.setObservations(
        jpype.JArray(jpype.JDouble)(source),
        jpype.JArray(jpype.JDouble)(destination),
    )

    assert (
        ours.computeAverageLocalOfObservations()
        == theirs.computeAverageLocalOfObservations()
    )


def test_two_dimensional_observations_match_a_hand_built_java_array(data):
    """2D arrays are the case JPype will not convert on its own."""
    source = data.normal(size=(1000, 2))
    destination = source + data.normal(size=(1000, 2)) * 0.5

    ours = jidt.MutualInfoCalculatorMultiVariateKraskov1()
    ours.initialise(2, 2)
    ours.setObservations(source, destination)

    theirs = jidt.MutualInfoCalculatorMultiVariateKraskov1()
    theirs.initialise(2, 2)
    theirs.setObservations(
        jpype.JArray(jpype.JDouble, 2)(source),
        jpype.JArray(jpype.JDouble, 2)(destination),
    )

    assert (
        ours.computeAverageLocalOfObservations()
        == theirs.computeAverageLocalOfObservations()
    )


def test_java_arrays_may_still_be_passed_by_the_user(data):
    """Code written against JIDT elsewhere must keep working unchanged."""
    x = data.integers(0, 2, 500).astype("int32")
    y = data.integers(0, 2, 500).astype("int32")

    calc = jidt.MutualInformationCalculatorDiscrete(2)
    calc.initialise()
    calc.addObservations(x, y)
    from_numpy = calc.computeAverageLocalOfObservations()

    calc.initialise()
    calc.addObservations(
        jpype.JArray(jpype.JInt)(x), jpype.JArray(jpype.JInt)(y)
    )

    assert calc.computeAverageLocalOfObservations() == from_numpy


def test_local_values_come_back_as_numpy(data):
    """A returned Java array becomes a numpy array of the same values.

    The values are those the toolkit gives when it is called directly
    through JPype, with this module out of the way.
    """
    source = data.normal(size=1000)
    destination = source + data.normal(size=1000) * 0.5

    calc = jidt.MutualInfoCalculatorMultiVariateKraskov1()
    calc.initialise(1, 1)
    calc.setObservations(source, destination)
    ours = calc.computeLocalOfPreviousObservations()

    java_class = jpype.JClass(
        "infodynamics.measures.continuous.kraskov."
        "MutualInfoCalculatorMultiVariateKraskov1"
    )
    unwrapped = java_class()
    unwrapped.initialise(1, 1)
    unwrapped.setObservations(
        jpype.JArray(jpype.JDouble)(source),
        jpype.JArray(jpype.JDouble)(destination),
    )
    theirs = unwrapped.computeLocalOfPreviousObservations()

    assert isinstance(ours, np.ndarray)
    assert not isinstance(theirs, np.ndarray)  # what JPype returns alone
    assert ours.shape == (1000,)
    assert np.array_equal(ours, np.asarray(theirs))


def test_objects_returned_from_java_stay_usable(data):
    """A returned calculator object keeps its fields and methods."""
    x = data.integers(0, 2, 1000)

    calc = jidt.MutualInformationCalculatorDiscrete(2)
    calc.initialise()
    calc.addObservations(x, data.integers(0, 2, 1000))
    calc.computeAverageLocalOfObservations()
    measurement = calc.computeSignificance(20)

    assert isinstance(measurement.pValue, float)
    assert isinstance(measurement.getMeanOfDistribution(), float)


def test_the_gui_can_be_launched_from_the_bundled_runtime():
    """Opening a window needs a display, so check what the launch needs.

    The GUI is started as "java -jar infodynamics.jar", which requires a
    java executable in the runtime and a main class in the jar's manifest.
    """
    assert (jidt.java_home() / "bin" / "java").exists()

    with zipfile.ZipFile(jidt.JAR) as jar:
        manifest = jar.read("META-INF/MANIFEST.MF").decode()

    assert "infodynamics.demos.autoanalysis.AutoAnalyserLauncher" in manifest


def test_matrix_utils_static_methods_take_python_data():
    utils = jidt.jclass("infodynamics.utils.MatrixUtils")
    assert utils.max([1.0, 7.0, 3.0]) == pytest.approx(7.0)
