# JIDT for Python

A lean Python wrapper around the [Java Information Dynamics
Toolkit](https://github.com/jlizier/jidt). Every calculator in the toolkit is
available under its own name, taking lists or numpy arrays. The
`infodynamics.jar` and a Java runtime come with the package.

```bash
pip3 install "git+https://github.com/jlizier/jidt.git#subdirectory=python"
```

The module lives in the `python` subdirectory of the toolkit's repository,
which is what `#subdirectory=python` names; `pip` clones the repository and
builds from there, so `git` is the only thing you need beforehand. Once the
package is published to PyPI, `pip3 install jidt` will do the same thing.

## Usage

Every JIDT class is available under its own name, directly on the module:

```python
import numpy as np
import jidt

calc = jidt.MutualInformationCalculatorDiscrete(2)
calc.initialise()
calc.addObservations(x, y)                     # lists or numpy arrays
print(calc.computeAverageLocalOfObservations())
```

```python
calc = jidt.MutualInfoCalculatorMultiVariateKraskov1()
calc.initialise(2, 2)                          # 2 source, 2 destination variables
calc.setObservations(source, destination)      # (samples, variables) numpy arrays
mi = calc.computeAverageLocalOfObservations()
locals_ = calc.computeLocalOfPreviousObservations()   # numpy array back
```

Through JPype alone, the same two lines require a JVM started with the jar on
its class path, and the calculator looked up in the Java package it lives in:

```python
jpype.startJVM(classpath=["/path/to/infodynamics.jar"])
package = jpype.JPackage("infodynamics.measures.discrete")
calc = package.MutualInformationCalculatorDiscrete(2)
```

The module does the startup and the lookup, and converts the data either way.

Names are unchanged from the Java toolkit, so any JIDT example, tutorial or
AutoAnalyser-generated code carries over. `dir(jidt)` lists everything
available; `jidt.jclass("infodynamics.utils.MatrixUtils")` reaches a class by
its full Java name.

## The AutoAnalyser GUI

```bash
python3 -m jidt
```

or from inside a session, where it opens and hands control straight back:

```python
import jidt

jidt.gui()
```

This is the same point-and-click application as the jar's, launched from the
runtime that came with the package. It runs in a JVM of its own, so `gui()`
returns immediately, the windows are unaffected by what you do in your session
and outlive it, and closing them leaves your session running. `gui()` hands
back the process, should you want to wait for the GUI to close or to stop it;
`python3 -m jidt` waits for it, which is what you want from a shell.

## What it does for you

- **Java runtime**: uses the one bundled by [`jdk4py`](https://pypi.org/project/jdk4py/).
  Set `JIDT_JAVA_HOME` (or `JAVA_HOME`) to use your own instead.
- **JVM startup**: automatic on first use, with `infodynamics.jar` on the
  classpath. Call `jidt.start("-Xmx4g")` first if you want your own JVM options.
- **64-bit check**: a 32-bit Python fails with a clear message rather than an
  obscure JVM error.
- **Data conversion**: lists and numpy arrays go in (2-D `(samples, variables)`
  arrays included); Java arrays come back as numpy arrays.

## Development

The jar is put into the package at build time: from `ant jar` in this checkout
if a JDK is available, otherwise from the full distribution attached to this
version's release on GitHub. The version itself comes from `build.xml`.

Install the module from this directory and run its tests with:

```bash
pip3 install -e ".[dev]"
pytest
```

These cover the wrapping only: that the jar and the runtime are found, that
the classes are reachable by name, and that data converts to and from Java
unchanged. The measures themselves are tested by the toolkit's own Java
suite, which ant runs as `ant junit`.
