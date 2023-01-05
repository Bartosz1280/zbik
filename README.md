# zbik
----
![Python](https://img.shields.io/badge/Python-3.10-blue)
![numpy](https://img.shields.io/badge/numpy-1.23-orange)
---
Zbik is an early-stage release of the project I did for python classes. The initial idea was to create a terminal-based program for sequence alignment by solving Rosalind's challenge.
Some bugs may appear. However, during the test, it did not show any errors.

The script requires at least `Python 3.10` to work. Regarding external libraries it only requires `NumPy` (one of the ideas was not to use `biopython` in the script).
 I prepared a file you can use to have a quick overview of the code.

`$ python run_demo.py`

`run_demo.py` uses files from subdirectories. I recommend `cd` to
a directory with the file and running it from there. In case, you run it from
some IDE, make sure about the correct `cwd`. I didn't have an opportunity to test it on an other OS than Linux. It should work the same way on other POSIX-compliant OS, but I am not sure if the way of path coding used in the script is compatible with Windows (i.e. `./sequences/`)

I didn't manage to implement sequences alignment before the deadline,
but I created a framework with an implementation that can solve four problems from Rosalind and has the potential for further development.

---
**TO-DO LIST:**

* Add remaining functionality to make Zbik a fully functional CLI tool
* Add distance matrices functionality
* Code refraction
